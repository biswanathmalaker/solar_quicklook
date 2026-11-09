#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
Common field of view: when did two or more instruments observe the same part of the Sun at the
same time? From all the cached metadata (no network), separate from the keyword search
(search_meta.py). Used by browse_web.py (the "common FOV" row and tab) and from the command line.

A result is a set of observations, one per instrument, going on at the same time whose fields of
view overlap: the common time, the common area as a share of the smallest field of view, and its
centre (Stonyhurst lat, lon). Observations:
  HRIEUV, HRI-Lya, PHI-HRT, SPICE   observing blocks (as the browser), FOV of the block's pointing
  FSI                               runs (gaps > 1 h), the 3.8 deg square around the Sun
  IRIS                              observations, the area the slit-jaw images cover (includes the raster)
  DKIST                             datasets, the bounding box
FOVs from different vantage points (Solar Orbiter, Earth) are compared through their footprints on
the surface, seen from above the footprints (so neither is hidden): a part off a limb counts at
that limb, so overlap in the corona between the two views is not counted. Same vantage point:
compared directly (off-limb parts included). Between Solar Orbiter and Earth only surface both see
at mu >= 0.1 counts (--mu; mu = cos of the viewing angle, 0 at the limb), so when Solar Orbiter is
round the far side a strip both see edge-on at a limb is not a common FOV; each result gives the mu
of its centre from each vantage point. IRIS full-disk mosaics (HCR gives the whole disk for hours,
IRIS sees one tile at a time) are left out unless --mosaics. The common area is measured on a
64 x 64 grid.

--without INST ...: instruments that must NOT be observing the common FOV. The times when one of
them shares any part of it (>= 1 % of the smallest FOV, same rules as above) are cut out of each
result; the pieces left, at least 5 min long (--left), are the results. One instrument is then
enough (e.g. HRIEUV blocks with no IRIS on them). The metadata knows the FOVs, not where a
feature is: an IRIS box inside the HRIEUV FOV but off the plume still excludes that time.

The index (time, observer, FOV corners per observation) is built per month from cache/soar,
cache/iris, cache/dkist into cache/fov_search/ (~3 min the first time), and a month is rebuilt when
its cached metadata is newer.

Run:  ./fov_search.py PHI-HRT IRIS --from 2025-10 --to 2025-11
      ./fov_search.py HRIEUV PHI-HRT SPICE --from 2025-10-01 --to 2025-10-31 --min 50 --max 100
      ./fov_search.py HRIEUV IRIS --mu 0.3 --mosaics
      ./fov_search.py HRIEUV --without IRIS --from 2024-03 --to 2024-05      # HRIEUV with no IRIS on its FOV
      ./fov_search.py HRIEUV PHI-HRT --without IRIS DKIST
      ./fov_search.py --build          # (re)build the index only
"""
import argparse
import shutil
import threading
from pathlib import Path

import numpy as np
import pandas as pd

import search_meta as sx
import soar_meta as sm

OUT_DIR = Path(__file__).resolve().parent
INDEX = OUT_DIR / "cache" / "fov_search"
GEOM = sx.GEOM
COLUMNS = ["instrument", "start", "end", "what"] + GEOM
_lock = threading.Lock()
_loaded = {}


def load(status=print):
    """All observations with their observer and FOV corners (one DataFrame, sorted by start)."""
    frames = []
    now = pd.Timestamp.now()
    for s, (_, _, first) in sx.SOURCES.items():
        for m in sm.months(first, now):
            df = sx.month_frame(s, m, status, index=INDEX, columns=COLUMNS, label="FOV index")
            if df is not None and len(df):
                frames.append(df)
    sig = tuple(id(f) for f in frames)
    with _lock:
        if _loaded.get("sig") == sig:
            return _loaded["df"]
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=COLUMNS)
    for c in GEOM:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=GEOM).sort_values("start", kind="stable").reset_index(drop=True)
    with _lock:
        _loaded.update(sig=sig, df=df)
    return df


GRID = 64               # the common area is measured on GRID x GRID points over the overlap of the boxes
MIN_PCT = 1.0           # a common area below 1 % of the smallest FOV is edges touching, not a common FOV
MIN_MU = 0.1            # default: the common area must be seen at least this face-on from both vantage points
MOSAIC = 1000.0         # arcsec: an IRIS FOV this wide is a full-disk mosaic (one tile at a time; HCR gives the whole disk)
EDGE_PTS = 12           # points per FOV edge (its footprint on the surface follows the curvature)
MIN_LEFT = pd.Timedelta("5min")   # "without": pieces of a result left between excluded observations shorter than this are dropped


def _rsun(r_au):
    """The Sun's apparent radius (arcsec) from r_au."""
    return np.degrees(np.arcsin(sm.RSUN_KM / (r_au * sm.AU_KM))) * 3600


PLATFORM = {"HRIEUV": "Solar Orbiter", "HRI-Lya": "Solar Orbiter", "PHI-HRT": "Solar Orbiter", "SPICE": "Solar Orbiter",
            "FSI": "Solar Orbiter", "IRIS": "Earth", "DKIST": "Earth"}


def _to_radius(x, y, mask, r):
    """Points in mask moved along their radius to r (arcsec from the disk centre)."""
    rho = np.hypot(x, y)
    f = np.where(mask & (rho > 0), r / np.where(rho > 0, rho, 1), 1.0)
    return x * f, y * f


def outline(r):
    """The FOV outline of an observation (its corners, each edge as EDGE_PTS points), arcsec in its own view."""
    c = np.array([[r[f"c{i}x"], r[f"c{i}y"]] for i in range(4)], float)
    t = np.linspace(0, 1, EDGE_PTS, endpoint=False)[:, None]
    return np.vstack([c[i] + t * (c[(i + 1) % 4] - c[i]) for i in range(4)])


def _view(r):
    return (r["obs_lat"], r["obs_lon"], r["obs_r"])


def _unit(lon, lat):
    lo, la = np.radians(lon), np.radians(lat)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def mu(lon, lat, view):
    """cos of the angle between the surface normal at Stonyhurst lon, lat and the direction to the
    observer view = (lat, lon, r_au) (exact for a finite distance; < 0 on the far side)."""
    n = _unit(np.atleast_1d(lon), np.atleast_1d(lat))
    o = _unit(np.atleast_1d(view[1]), np.atleast_1d(view[0]))[0] * view[2] * sm.AU_KM / sm.RSUN_KM
    d = o - n
    return np.sum(n * d, axis=1) / np.linalg.norm(d, axis=1)


def footprint(r, min_mu=0.0):
    """The FOV outline of an observation on the solar surface: unit vectors (n, 3), Stonyhurst.
    Only the part its observer sees at mu >= min_mu: points beyond are taken at that circle (the
    part off the limb at the limb when min_mu = 0). Also its centre (unit vector): the mean
    direction, or for a FOV around the whole disk (FSI) the observer's sub-point."""
    pts, src = outline(r), _view(r)
    rs = _rsun(src[2]) * min(0.999, np.sqrt(1 - min_mu ** 2))
    x, y = _to_radius(pts[:, 0], pts[:, 1], np.hypot(pts[:, 0], pts[:, 1]) > rs, rs)
    lon, lat = sm.hpc_to_hgs(x, y, *src)
    ok = np.isfinite(lon) & np.isfinite(lat)
    v = _unit(lon[ok], lat[ok])
    m = v.mean(axis=0) if len(v) else np.zeros(3)
    centre = m / np.linalg.norm(m) if np.linalg.norm(m) > 0.5 else _unit(src[1], src[0])[0]
    return v, centre


def overlap(rows, cache, min_mu=MIN_MU):
    """The part of the Sun common to the FOVs of these observations, or None: (% of the smallest
    FOV, lat, lon of its centre, index of the smallest). Same platform (Solar Orbiter; or Earth:
    IRIS, DKIST): compared in that view (off-limb parts included). Else: their footprints on the surface, compared in a view
    from above the footprints (so neither is hidden; off-limb parts count at the limb), only the
    parts each vantage point sees at mu >= min_mu."""
    views = [_view(r) for r in rows]
    if len({PLATFORM[r["instrument"]] for r in rows}) == 1:
        # one spacecraft (or Earth) at the common time; each block's pointing and position are its
        # medians, so outlines are taken in units of the Sun's apparent radius from that position
        polys = [outline(r) / _rsun(v[2]) for r, v in zip(rows, views)]
        res = common_area(polys)
        if res is None:
            return None
        area, cx, cy = res
        rs = _rsun(views[0][2])
        lon, lat = sm.hpc_to_hgs(np.array([cx * rs]), np.array([cy * rs]), *views[0])
        lon, lat = float(lon[0]), float(lat[0])
    else:
        feet = []
        for r in rows:
            k = (id(r), min_mu)
            if k not in cache:
                cache[k] = footprint(r, min_mu)
            feet.append(cache[k])
        if any(len(v) < 3 for v, _ in feet):
            return None
        c = np.sum([ctr for _, ctr in feet], axis=0)
        if np.linalg.norm(c) < 1e-6:
            return None
        c = c / np.linalg.norm(c)
        e1 = np.cross([0.0, 0.0, 1.0], c)
        e1 = e1 / np.linalg.norm(e1) if np.linalg.norm(e1) > 1e-9 else np.array([1.0, 0.0, 0.0])
        e2 = np.cross(c, e1)
        polys = []
        for v, _ in feet:   # orthographic view along c; points beyond its horizon put on the horizon
            x, y = _to_radius(v @ e1, v @ e2, v @ c < 0, 1.0)
            polys.append(np.column_stack([x, y]))
        res = common_area(polys)
        if res is None:
            return None
        area, cx, cy = res
        q = cx * e1 + cy * e2 + np.sqrt(max(0.0, 1 - cx * cx - cy * cy)) * c
        lon, lat = float(np.degrees(np.arctan2(q[1], q[0]))), float(np.degrees(np.arcsin(np.clip(q[2], -1, 1))))
    sizes = [_area(p) for p in polys]
    k = int(np.argmin(sizes))
    pct = min(100.0, 100 * area / sizes[k]) if sizes[k] > 0 else np.nan
    if not pct >= MIN_PCT:      # edges that only touch
        return None
    return pct, lat, lon, k


def _area(p):
    x, y = p[:, 0], p[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def common_area(polys):
    """The area (arcsec^2) and centre (x, y) common to all the outlines (same view), or None."""
    from matplotlib.path import Path as MPath
    lo, hi = np.max([p.min(axis=0) for p in polys], axis=0), np.min([p.max(axis=0) for p in polys], axis=0)
    if np.any(lo >= hi):
        return None
    step = (hi - lo) / GRID
    gx, gy = np.meshgrid(lo[0] + step[0] * (np.arange(GRID) + 0.5), lo[1] + step[1] * (np.arange(GRID) + 0.5))
    pts = np.column_stack([gx.ravel(), gy.ravel()])
    inside = np.ones(len(pts), bool)
    for p in polys:
        inside &= MPath(p).contains_points(pts)
        if not inside.any():
            return None
    c = pts[inside].mean(axis=0)
    return inside.sum() * step[0] * step[1], c[0], c[1]


def _excluded(members, s, e, xrows, cache, min_mu):
    """The times in [s, e] when an excluded observation (xrows: (starts, ends, rows) per instrument)
    shares part of the members' common FOV: a list of (start, end, instrument), merged per instrument."""
    out = []
    for inst, (starts, ends, rows, maxdur) in xrows.items():
        lo = np.searchsorted(starts, np.datetime64(s - maxdur), "left")
        hi = np.searchsorted(starts, np.datetime64(e), "right")
        for j in range(lo, hi):
            if ends[j] < np.datetime64(s):
                continue
            if overlap(members + [rows[j]], cache, min_mu) is not None:
                out.append((max(s, rows[j]["start"]), min(e, rows[j]["end"]), inst))
    return sorted(out)


def _left(s, e, cut, min_left):
    """[s, e] minus the cut intervals: the pieces at least min_left long."""
    pieces, t = [], s
    for a, b, _ in sorted(cut):
        if a > t:
            pieces.append((t, a))
        t = max(t, b)
    if e > t:
        pieces.append((t, e))
    return [(a, b) for a, b in pieces if b - a >= min_left]


def together_fov(hits, instruments, status=print, min_mu=MIN_MU, exclude=None, min_left=MIN_LEFT, without=None):
    """Like together(), but each result is a set of observations, one per instrument, going on
    at the same time AND overlapping on the Sun: instrument, start, end (common time), what (the
    common area as a share of the smallest FOV, its centre, and each observation), overlap (%),
    lat, lon (centre of the common area, Stonyhurst; NaN off the limb). See overlap().
    exclude: observations of other instruments (same columns) that must NOT be on the common FOV:
    the times when one of them shares part of it are cut out of the result, and the pieces left at
    least min_left long are the results (one instrument alone is then allowed). without: the
    excluded instruments' names (listed in each result even if they have no observations)."""
    instruments = [i for i in sx.INSTRUMENTS if i in instruments]
    have = hits.dropna(subset=GEOM)
    groups = [have[have.instrument == i].sort_values("start").reset_index(drop=True) for i in instruments]
    if any(not len(g) for g in groups):
        return pd.DataFrame(columns=["instrument", "start", "end", "what", "overlap", "lat", "lon"])
    # partial sets: (row numbers, start, end); each new member is checked against all the earlier ones
    rows = [g.to_dict("records") for g in groups]
    cache = {}
    partial = [([k], r["start"], r["end"]) for k, r in enumerate(rows[0])]
    for n, g in enumerate(groups[1:], 1):
        status(f"common FOV: {instruments[0]} .. {instruments[n]}, {len(partial)} candidate(s) ...")
        starts, ends = g.start.values, g.end.values
        maxdur = (g.end - g.start).max()
        new = []
        for ids, s, e in partial:
            lo = np.searchsorted(starts, np.datetime64(s - maxdur), "left")
            hi = np.searchsorted(starts, np.datetime64(e), "right")
            for j in range(lo, hi):
                if ends[j] < np.datetime64(s):
                    continue
                r = rows[n][j]
                if overlap([rows[m][i] for m, i in enumerate(ids)] + [r], cache, min_mu) is None:
                    continue
                new.append((ids + [j], max(s, r["start"]), min(e, r["end"])))
        partial = new
    xrows, names = {}, []
    if exclude is not None:
        exclude = exclude.dropna(subset=GEOM)
        names = [i for i in sx.INSTRUMENTS if i in set(exclude.instrument) | set(without or ())]
        for inst in names:
            g = exclude[exclude.instrument == inst].sort_values("start").reset_index(drop=True)
            if len(g):
                xrows[inst] = (g.start.values, g.end.values, g.to_dict("records"), (g.end - g.start).max())
        status(f"common FOV: {len(partial)} set(s), cutting out the times with {', '.join(names)} on them ...")
    out = []
    for ids, s, e in partial:
        members = [rows[m][i] for m, i in enumerate(ids)]
        pct, lat, lon, k = overlap(members, cache, min_mu)
        where = f"centre lat {lat:+.0f}°, lon {lon:+.0f}°" if np.isfinite(lat) else "centre off the limb"
        if np.isfinite(lat):   # how face-on each vantage point sees the centre
            seen = {}
            for x in members:
                seen.setdefault(PLATFORM[x["instrument"]], float(mu(lon, lat, _view(x))[0]))
            where += ", seen at " + ", ".join(f"mu {m:.2f} from {pl}" for pl, m in seen.items())
        if len(instruments) == 1:
            lines = [f"{instruments[0]} FOV {where}"]
        else:
            lines = [f"common FOV: {pct:.0f}% of the {instruments[k]} FOV, {where}"]
        lines += [f"{inst}: " + " ".join(rows[m][i]["what"].split()) for m, (inst, i) in enumerate(zip(instruments, ids))]
        pieces = [(s, e)]
        if names:
            cut = _excluded(members, s, e, xrows, cache, min_mu)
            pieces = _left(s, e, cut, min_left)
            note = f"without {', '.join(names)} on this FOV"
            if cut:
                note += "; left out, when on it: " + ", ".join(
                    f"{inst} {a:%H:%M}-{b:%H:%M}" for a, b, inst in cut[:6]) + (" ..." if len(cut) > 6 else "")
            lines.append(note)
        for a, b in pieces:
            out.append(dict(instrument=" + ".join(instruments), start=a, end=b, what="\n".join(lines),
                            overlap=pct, lat=lat, lon=lon))
    df = pd.DataFrame(out, columns=["instrument", "start", "end", "what", "overlap", "lat", "lon"])
    if len(df):   # several datasets of one area at once (e.g. DKIST, one per line): one result, all of them listed
        key = [df.start.dt.floor("min").rename("_s"), df.end.dt.floor("min").rename("_e"),
               (df.overlap / 5).round(0).rename("_o"), df.lat.round(0).rename("_la"), df.lon.round(0).rename("_lo")]
        df = (df.groupby(key, sort=False, dropna=False)
                .agg(instrument=("instrument", "first"), start=("start", "min"), end=("end", "max"),
                     what=("what", _merge), overlap=("overlap", "first"), lat=("lat", "first"), lon=("lon", "first"))
                .reset_index(drop=True))
    return df[["instrument", "start", "end", "what", "overlap", "lat", "lon"]].sort_values("start").reset_index(drop=True)


def _merge(whats):
    """The text of equivalent results: each distinct line once, in order."""
    return "\n".join(dict.fromkeys(line for w in whats for line in w.split("\n")))


def is_mosaic(df):
    """IRIS full-disk mosaics: HCR gives the whole disk, but IRIS sees one tile at a time."""
    w = df[["c0x", "c1x", "c2x", "c3x"]].max(axis=1) - df[["c0x", "c1x", "c2x", "c3x"]].min(axis=1)
    h = df[["c0y", "c1y", "c2y", "c3y"]].max(axis=1) - df[["c0y", "c1y", "c2y", "c3y"]].min(axis=1)
    return (df.instrument == "IRIS") & ((w > MOSAIC) | (h > MOSAIC))


def find(instruments, t0=None, t1=None, min_overlap=0.0, status=print, min_mu=MIN_MU, mosaics=False,
         without=(), min_left=MIN_LEFT):
    """Sets of observations of all the instruments (two or more; one is enough with without) at
    the same time on a common part of the Sun, in [t0, t1], whose common area is at least
    min_overlap % of the smallest FOV and, between Solar Orbiter and Earth, seen at mu >= min_mu
    from both. IRIS full-disk mosaics only with mosaics=True. without: instruments that must not
    be observing any part of that common area (the times they are, are cut out)."""
    without = [i for i in without if i not in instruments]
    df = load(status)
    df = df[df.instrument.isin(list(instruments) + without)]
    if not mosaics:
        df = df[~is_mosaic(df)]
    if t0 is not None:
        df = df[df.end >= pd.Timestamp(t0)]
    if t1 is not None:
        df = df[df.start <= pd.Timestamp(t1)]
    out = together_fov(df[df.instrument.isin(instruments)], instruments, status, min_mu,
                       exclude=df[df.instrument.isin(without)] if without else None, min_left=min_left,
                       without=without)
    return out[out.overlap >= min_overlap].reset_index(drop=True) if min_overlap else out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inst", nargs="*", metavar="INST", help="two or more of " + " ".join(sx.INSTRUMENTS))
    ap.add_argument("--from", dest="t0", help="observations ending after this (UTC, e.g. 2025-10)")
    ap.add_argument("--to", dest="t1", help="observations starting before this")
    ap.add_argument("--min", type=float, default=0, help="common area at least this %% of the smallest FOV (default 0)")
    ap.add_argument("--mu", type=float, default=MIN_MU, help="Solar Orbiter vs Earth: common area seen at least this "
                    f"face-on (mu) from both (default {MIN_MU}; polar regions are ~0.2-0.35 from Earth)")
    ap.add_argument("--mosaics", action="store_true", help="include IRIS full-disk mosaics (whole disk for hours, "
                    "one tile at a time)")
    ap.add_argument("--without", nargs="+", default=[], metavar="INST", help="instruments that must NOT observe the "
                    "common FOV: the times they do are cut out (with this, one instrument is enough)")
    ap.add_argument("--left", type=float, default=MIN_LEFT.total_seconds() / 60, help="with --without: pieces "
                    f"left shorter than this many minutes are dropped (default {MIN_LEFT.total_seconds() / 60:g})")
    ap.add_argument("--max", type=int, default=50, help="results printed (default 50; 0: all)")
    ap.add_argument("--build", action="store_true", help="(re)build the index from the cache only")
    args = ap.parse_args()
    bad = [i for i in args.inst + args.without if i not in sx.INSTRUMENTS]
    if bad:
        ap.error(f"unknown instrument(s) {' '.join(bad)}; choose from {' '.join(sx.INSTRUMENTS)}")
    both = [i for i in args.without if i in args.inst]
    if both:
        ap.error(f"{' '.join(both)} both wanted and --without")
    if args.build or len(args.inst) < (1 if args.without else 2):
        df = load()
        print(f"FOV index: {len(df)} observations with a FOV, " + ", ".join(f"{k} {v}" for k, v in df.instrument.value_counts().items()))
        if not args.build:
            print("give two or more instruments to search (one is enough with --without)")
        return
    hits = find(args.inst, args.t0, args.t1, args.min, status=lambda s: None, min_mu=args.mu, mosaics=args.mosaics,
                without=args.without, min_left=pd.Timedelta(minutes=args.left))
    print(f"{len(hits)} sets of {' + '.join(args.inst)} observations at the same time on a common FOV"
          + (f", without {' '.join(args.without)} on it" if args.without else "")
          + (f" (first {args.max})" if 0 < args.max < len(hits) else ""))
    width = shutil.get_terminal_size((160, 20)).columns
    for _, r in (hits.head(args.max) if args.max else hits).iterrows():
        line = f"{r.start:%Y-%m-%d %H:%M} .. {r.end:%Y-%m-%d %H:%M}  " + r.what.replace("\n", "  ##  ")
        print(line[:width - 1])


if __name__ == "__main__":
    main()
