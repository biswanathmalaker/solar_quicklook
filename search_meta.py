#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
Keyword search over all the cached metadata (no network): Solar Orbiter observing blocks and FSI
runs (SOAR), IRIS observations (LMSAL HCR) and DKIST datasets. Used by browse_web.py (the
"search" tab) and from the command line.

Matching: case-insensitive; keywords separated by commas, any of them (OR); a keyword of
several words is a phrase (those words in that order). Words match at a word start ("hole" finds
"holes", not "whole"); "-", "_", "/" and punctuation count as spaces, so "polar observations"
finds the SOOP R_SMALL_MRES_MCAD_Polar-Observations. With two or more instruments, the result
is the times when all of them have a matching observation going on at once (AND): one row per
common interval, with the overlapping observations of each. (Observed together on a common
part of the Sun, without keywords: fov_search.py.)

What is searched:
  HRIEUV, HRI-Lya,  instrument, product, SOOP name, the SOOP summary's description, coordinators
  PHI-HRT, SPICE    and notes (cache/soops), the SOOP's ESA page when cached, SPICE study names,
                    file descriptors; one result per observing block
  FSI               "FSI 174" and the SOOP names of the run (runs split at gaps > 1 h)
  IRIS              OBSID, aim, planners, programme, target, description, spectral lines, slit-jaw channels
  DKIST             instrument, proposal, experiment and its description, target types, lines, dataset ID

The index is built per month from cache/soar, cache/iris and cache/dkist into cache/search/ (the
first time ~3 min, mostly the Solar Orbiter blocks), and a month is rebuilt when its cached
metadata is newer. Months not in the cache are not searched (fetch them with --all first).

Run:  ./search_meta.py "coronal hole, plume"
      ./search_meta.py "coronal hole, polar" --inst IRIS SPICE --from 2023-01 --to 2024-12 --max 100
      ./search_meta.py --build          # (re)build the index only
"""
import argparse
import json
import re
import shutil
import threading
from pathlib import Path

import numpy as np
import pandas as pd

import dkist_meta as dm
import iris_meta as im
import soar_meta as sm
import soops

OUT_DIR = Path(__file__).resolve().parent
INDEX = OUT_DIR / "cache" / "search" / "v1"      # bump the version when the index format changes
INSTRUMENTS = ["HRIEUV", "HRI-Lya", "PHI-HRT", "SPICE", "FSI", "IRIS", "DKIST"]
LABELS = {"hri174": "HRIEUV", "hriopn": "HRIEUV", "hrilya": "HRI-Lya", "phi-hrt": "PHI-HRT", "spice": "SPICE"}
CORNERS = ["c0x", "c0y", "c1x", "c1y", "c2x", "c2y", "c3x", "c3y"]
GEOM = ["obs_lat", "obs_lon", "obs_r"] + CORNERS      # observer (Stonyhurst deg, au); FOV corners (arcsec, its view)
COLUMNS = ["instrument", "start", "end", "soop", "what", "text"]   # the keyword index
ROW_COLUMNS = COLUMNS + GEOM        # what the month builders give (the geometry: fov_search.py)
FSI_HALF = 6840.0       # arcsec: half of FSI's 3.8 deg field of view
FSI_RUN_GAP = pd.Timedelta("1h")
DKIST_TARGETS = {"coronalhole": "coronal hole", "activecorona": "active corona", "quietcorona": "quiet corona",
                 "quietsun": "quiet sun"}   # DKIST targetTypes are joined words
_lock = threading.Lock()
_frames = {}            # index file -> (mtime, DataFrame)
_loaded = {}            # sources -> (signature, DataFrame)


def norm(text):
    """Lower case, punctuation and -_/ as single spaces (index text and query words alike)."""
    return " ".join(re.sub(r"[\W_]+", " ", str(text).lower()).split())


def fmt_s(s):
    if not np.isfinite(s):
        return "?"
    return f"{s:.0f} s" if s < 120 else f"{s / 60:.0f} min" if s < 7200 else f"{s / 3600:.1f} h"


def _clean(v):
    return "" if v is None or (isinstance(v, float) and np.isnan(v)) or str(v).lower() == "nan" else str(v)


# ---------------------------------------------------------------- one month, per source

def _end(m):
    return m + pd.offsets.MonthBegin(1) - pd.Timedelta("1us")


def earth_observer(times):
    """Earth's Stonyhurst lat (B0, deg), lon (0) and distance (au) at these times."""
    from astropy.time import Time
    from sunpy.coordinates import sun
    if not len(times):
        return np.zeros(0), np.zeros(0), np.zeros(0)
    t = Time(list(pd.to_datetime(pd.Series(times)).dt.tz_localize(None).dt.to_pydatetime()))
    return sun.B0(t).deg, np.zeros(len(t)), sun.earth_distance(t).to("au").value


def box(x0, y0, x1, y1):
    return dict(c0x=x0, c0y=y0, c1x=x1, c1y=y0, c2x=x1, c2y=y1, c3x=x0, c3y=y1)


def solo_sources(m):
    """The cached SOAR files of month m, or None if the month is not (fully) cached."""
    need = [sm.CACHE / f"items_{p}_{m:%Y-%m}.csv.gz" for p in list(sm.POINTED) + ["fsi174"]]
    if not all(p.is_file() for p in need):
        return None
    return sorted(sm.CACHE.glob(f"*_{m:%Y-%m}.csv.gz"))


def solo_month(m):
    """Observing blocks (as browse.MetaStore) and FSI 174 runs of month m, from the cache."""
    rows = []
    f = sm.files(m, _end(m), sm.POINTED)
    b, f = sm.blocks(f) if len(f) else (pd.DataFrame(), f)
    if len(f) and (f["product"] == "spice").any():
        modes = sm.spice_modes(m)
        if len(modes):
            f = f.merge(modes, on="data_item_id", how="left")
    for _, r in (b.iterrows() if len(b) else []):
        label = LABELS.get(r["product"])
        if label is None:
            continue
        studies = []
        if "study" in f:
            studies = sorted(set(f.loc[f.block == r.block, "study"].dropna().astype(str)))
        n = int(r.n)
        what = [f"{n} files, one every {fmt_s(r.cadence_s)}" if n > 1 else "1 file"]
        if np.isfinite(r.get("lat_min", np.nan)):
            what.append(f"lat {r.lat_min:+.0f}..{r.lat_max:+.0f}°")
        what.append(f"SOOP {r.soop}")
        if studies:
            what.append("study " + ", ".join(studies))
        text = [label, r["product"], r.soop, r.descriptors] + studies
        rows.append(dict(instrument=label, start=r.start, end=r.end, soop=r.soop, what="; ".join(what),
                         text=norm(" ".join(map(str, text))), obs_lat=r.solo_lat, obs_lon=r.solo_lon, obs_r=r.solo_r,
                         **{c: r[c] for c in CORNERS}))
    fsi = sm.files(m, _end(m), ["fsi174"], headers=False)
    # Solar Orbiter's position for FSI (its files here have none): from the pointed files of the month
    pos = f.dropna(subset=["solo_lat", "solo_lon", "solo_r"]).sort_values("start") if len(f) and "solo_lat" in f else []
    def where(t):
        if not len(pos):
            return dict(obs_lat=np.nan, obs_lon=np.nan, obs_r=np.nan)
        x = pos.start.values.astype("datetime64[ns]").astype(np.int64)
        v = np.datetime64(t, "ns").astype(np.int64)
        return {k: float(np.interp(v, x, pos[c])) for k, c in (("obs_lat", "solo_lat"), ("obs_lon", "solo_lon"),
                                                                 ("obs_r", "solo_r"))}
    if len(fsi):
        fsi = fsi.sort_values("start")
        run = (fsi.start.diff() > FSI_RUN_GAP).cumsum()
        for _, x in fsi.groupby(run):
            cad = x.start.diff().dt.total_seconds().median()
            ss = sorted(set(x.soop) - {"none"})
            rows.append(dict(instrument="FSI", start=x.start.iloc[0], end=x.start.iloc[-1], soop=";".join(ss),
                             what=f"FSI 174: {len(x)} images, one every {fmt_s(cad)}" + (f"; SOOP {', '.join(ss)}" if ss else ""),
                             text=norm("FSI 174 " + " ".join(ss)), **where(x.start.iloc[0]),
                             **box(-FSI_HALF, -FSI_HALF, FSI_HALF, FSI_HALF)))
    return pd.DataFrame(rows, columns=ROW_COLUMNS)


def iris_sources(m):
    p = im.CACHE / f"hcr_{m:%Y-%m}.json.gz"
    return [p] if p.is_file() else None


def iris_month(m):
    """IRIS observations of month m (from the cached HCR answer; not kept in iris_meta.MONTHS)."""
    rows = []
    for e in im.query_month(m, verbose=False):
        r = im._row(e)
        aim, goal, planners = r["aim"], r["goal"], r["planners"]
        what = f"OBSID {r['obsid']}: {aim or '(no aim given)'} | {goal} | {planners} | ({r['xcen']:+.0f}\", {r['ycen']:+.0f}\")"
        text = [r["obsid"], aim, planners, goal, r["target"], r["description"], r["lines"], "sji " + r["sji_channels"]]
        c = im.corners(pd.Series(r), "sji")      # the area the slit-jaw images cover (includes the raster)
        if c is None:
            c = im.corners(pd.Series(r), "raster")
        geom = {k: v for k, v in zip(CORNERS, c.ravel())} if c is not None else {}
        rows.append(dict(instrument="IRIS", start=r["start"], end=r["end"], soop="", what=what,
                         text=norm(" ".join(map(_clean, ["IRIS"] + text))), **geom))
    df = pd.DataFrame(rows, columns=ROW_COLUMNS)
    df["obs_lat"], df["obs_lon"], df["obs_r"] = earth_observer(df.start)
    return df


def dkist_sources(m):
    p = dm.CACHE / f"datasets_{m:%Y-%m}.csv.gz"
    return [p] if p.is_file() else None


def dkist_month(m):
    d = dm.query_month(m, verbose=False)
    rows = []
    for _, r in (d.iterrows() if len(d) else []):
        g = lambda k: _clean(r.get(k))
        desc = g("experimentDescription")
        what = (f"{g('instrumentName')} {g('spectralLines')}; targets {g('targetTypes') or '?'}; "
                f"{g('primaryProposalId')} / {g('primaryExperimentId')}; {g('datasetId')}" + (f": {desc}" if desc else ""))
        targets = g("targetTypes")
        text = ["DKIST"] + [g(k) for k in ("instrumentName", "primaryProposalId", "primaryExperimentId",
                                           "experimentDescription", "spectralLines", "datasetId")]
        text += [targets] + [DKIST_TARGETS.get(t.strip(), "") for t in targets.split(";")]
        rows.append(dict(instrument="DKIST", start=pd.Timestamp(r.startTime), end=pd.Timestamp(r.endTime), soop="",
                         what=what, text=norm(" ".join(text)), **box(*dm.parse_box(r.boundingBox))))
    df = pd.DataFrame(rows, columns=ROW_COLUMNS)
    df["obs_lat"], df["obs_lon"], df["obs_r"] = earth_observer(df.start)
    return df


SOURCES = {"solo": (solo_sources, solo_month, sm.FIRST_MONTH),
           "iris": (iris_sources, iris_month, im.FIRST_MONTH),
           "dkist": (dkist_sources, dkist_month, dm.FIRST_MONTH)}


# ---------------------------------------------------------------- the index

def month_frame(source, m, status=print, index=None, columns=None, label="search index"):
    """The index of one source and month (cache/search/; or index, keeping these columns), rebuilt
    if its metadata is newer; None if the month's metadata is not in the cache."""
    index = INDEX if index is None else index
    columns = COLUMNS if columns is None else columns
    sources, build, _ = SOURCES[source]
    src = sources(m)
    if not src:
        return None
    path = index / f"{source}_{m:%Y-%m}.csv.gz"
    newest = max(p.stat().st_mtime for p in src)
    if not path.is_file() or path.stat().st_mtime < newest:
        status(f"{label}: {source} {m:%Y-%m} ...")
        df = build(m)[columns]
        index.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        df.to_csv(tmp, index=False, compression="gzip")
        tmp.replace(path)
    mtime = path.stat().st_mtime
    with _lock:
        if path in _frames and _frames[path][0] == mtime:
            return _frames[path][1]
    df = pd.read_csv(path, keep_default_na=False, dtype=str)
    df["start"] = pd.to_datetime(df.start, format="mixed")
    df["end"] = pd.to_datetime(df.end, format="mixed")
    with _lock:
        _frames[path] = (mtime, df)
    return df


def soop_texts():
    """SOOP full name -> its summary description, coordinators, notes and cached ESA page text (normalised)."""
    out = {}
    try:
        t = soops.table()
    except Exception:
        return out
    for name, g in t.groupby("soop_full"):
        parts = []
        for c in ("soop", "description", "coordinators", "notes"):
            parts += list(dict.fromkeys(_clean(v) for v in g[c]))
        page = soops.PAGES / f"{name}.json"
        if page.is_file():
            p = json.loads(page.read_text())
            parts.append(p.get("description", ""))
            parts += [" ".join(map(str, row)) for row in p.get("instruments", [])]
        out[str(name)] = norm(" ".join(parts))
    return out


def load(sources=tuple(SOURCES), status=print):
    """All the indexed observations (one DataFrame, sorted by start)."""
    frames = []
    now = pd.Timestamp.now()
    for s in sources:
        for m in sm.months(SOURCES[s][2], now):
            df = month_frame(s, m, status)
            if df is not None and len(df):
                frames.append(df)
    sig = (tuple(id(f) for f in frames), soops.CSV.stat().st_mtime if soops.CSV.is_file() else 0,
           len(list(soops.PAGES.glob("*.json"))) if soops.PAGES.is_dir() else 0)
    with _lock:
        if sources in _loaded and _loaded[sources][0] == sig:
            return _loaded[sources][1]
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=COLUMNS)
    texts = soop_texts()
    if texts and len(df):   # the SOOP descriptions, added here (they change apart from the SOAR months)
        extra = df.soop.map(lambda s: " ".join(texts.get(x, "") for x in s.split(";") if x)).fillna("")
        df["text"] = df.text + " " + extra
    df["text"] = " " + df.text          # words then start after a space (search)
    df = df.sort_values("start", kind="stable").reset_index(drop=True)
    with _lock:
        _loaded[sources] = (sig, df)
    return df


def terms(query):
    """The query's keywords: comma separated, each a word or a phrase (words in that order)."""
    return [t for t in (norm(x) for x in query.replace('"', " ").split(",")) if t]


def search(query, instruments=None, t0=None, t1=None, index=None, status=print):
    """Observations of the instruments whose text contains any of the keywords (OR):
    instrument, start, end, what."""
    df = load(status=status) if index is None else index
    if instruments:
        df = df[df.instrument.isin(instruments)]
    if t0 is not None:
        df = df[df.end >= pd.Timestamp(t0)]
    if t1 is not None:
        df = df[df.start <= pd.Timestamp(t1)]
    keys = [" " + t for t in terms(query)]      # at a word start: "hole" finds "holes", not "whole"
    if keys:            # no keywords: every observation
        df = df[[any(k in text for k in keys) for text in df.text]]
    return df[["instrument", "start", "end", "what"]].reset_index(drop=True)


def _union(df):
    """Sorted, merged [start, end] intervals covered by these observations."""
    out = []
    for s, e in sorted(zip(df.start, df.end)):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def _intersect(a, b):
    out, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s <= e:
            out.append([s, e])
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def together(hits, instruments, per_instrument=2):
    """The times when every one of the instruments has a matching observation going on (AND):
    one row per common interval: instrument ("IRIS + SPICE"), start, end (of the common
    interval), what (the overlapping observations of each instrument)."""
    instruments = [i for i in INSTRUMENTS if i in instruments]
    common = None
    for inst in instruments:
        u = _union(hits[hits.instrument == inst])
        common = u if common is None else _intersect(common, u)
    rows = []
    for s, e in common or []:
        parts = []
        for inst in instruments:
            g = hits[(hits.instrument == inst) & (hits.start <= e) & (hits.end >= s)]
            more = f" (+{len(g) - per_instrument} more)" if len(g) > per_instrument else ""
            parts.append(f"{inst}: " + " || ".join(" ".join(w.split()) for w in g.what.head(per_instrument)) + more)
        rows.append(dict(instrument=" + ".join(instruments), start=s, end=e, what="\n".join(parts)))
    return pd.DataFrame(rows, columns=["instrument", "start", "end", "what"])


def find(query, instruments, t0=None, t1=None, index=None, status=print):
    """The search as the browser does it: keywords OR-ed (commas); with two or more instruments
    only the times all of them observed together (together()), else their matching observations."""
    hits = search(query, instruments, t0, t1, index, status)
    return together(hits, instruments) if instruments and len(instruments) > 1 else hits


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("query", nargs="?", help='keywords separated by commas, any of them (OR), e.g. "coronal hole, plume"')
    ap.add_argument("--inst", nargs="+", choices=INSTRUMENTS, metavar="INST",
                    help="instruments: " + " ".join(INSTRUMENTS) + "; with two or more, only the times when all of "
                         "them have a matching observation at once (AND); default: all, each on its own")
    ap.add_argument("--from", dest="t0", help="observations ending after this (UTC, e.g. 2023-01)")
    ap.add_argument("--to", dest="t1", help="observations starting before this")
    ap.add_argument("--max", type=int, default=50, help="results printed (default 50; 0: all)")
    ap.add_argument("--build", action="store_true", help="(re)build the index from the cache only")
    args = ap.parse_args()
    if args.build or not args.query:
        df = load()
        print(f"search index: {len(df)} observations, " + ", ".join(f"{k} {v}" for k, v in df.instrument.value_counts().items()))
        return
    hits = find(args.query, args.inst, args.t0, args.t1, status=lambda s: None)
    kind = "times observed together" if args.inst and len(args.inst) > 1 else "observations"
    print(f"{len(hits)} {kind} match {' OR '.join(terms(args.query))}" + (f" (first {args.max})" if 0 < args.max < len(hits) else ""))
    width = shutil.get_terminal_size((160, 20)).columns
    for _, r in (hits.head(args.max) if args.max else hits).iterrows():
        line = f"{r.instrument:8s} {r.start:%Y-%m-%d %H:%M} .. {r.end:%Y-%m-%d %H:%M}  " + r.what.replace("\n", "  ##  ")
        print(line[:width - 1])


if __name__ == "__main__":
    main()
