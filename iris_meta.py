#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
IRIS observations in a time range, from the LMSAL Heliophysics Coverage Registry (HCR), the
database behind iris.lmsal.com/search (no login):

    https://www.lmsal.com/hek/hcr?cmd=search-events3&outputformat=json&instrument=IRIS&startTime=...&stopTime=...

One record per observation (OBSID run): the aim (sciObjectives, e.g. "C1: North Pole" or
"Solar Orbiter EUI/SPICE coordination"), planner, the programme ("goal", e.g. "Very large
coarse 64-step raster 126x175 64s"), pointing (xcen, ycen: arcsec, Earth view, at the
start), roll, FOVs, spectral windows, and per raster / slit-jaw (SJI) channel the cadence,
number of steps, step size, exposure, image count and the level-2 data links and sizes.
The answer of a month is kept as it came (cache/iris/hcr_YYYY-MM.json.gz): queried once,
again only with --refresh (e.g. for recent days whose level-2 entries were not ready yet).

FOVs (all centred on xcen, ycen, rotated by the roll):
  raster   raster_w x raster_h: the area the spectrograph slit covers (steps x step size,
           slit length); the archive's own footprint (hpc_geom_inner) is the same box
  SJI      cover_w x cover_h: the area the slit-jaw images cover over the raster (the SJI
           field moves with the slit, so it is about raster_w + sji_w wide); one SJI
           image is sji_w x sji_h

    import iris_meta as im
    d = im.observations("2025-10-10", "2025-10-17")      # one row per observation
    im.corners(d.iloc[0], "raster")                      # 4 x 2 arcsec

Run:  ./iris_meta.py 2025-10-10 2025-10-17               # table
      ./iris_meta.py 2025-10-10 2025-10-17 --pole        # only |ycen| > 800" (polar)
      ./iris_meta.py 2025-10-13T20:00 2025-10-14T00:00 --long    # with aim, planner, channels
      ./iris_meta.py --all       # fetch every month since the start (2013-07) not yet in the cache
                                 # (~40 s a month, ~0.5 MB; months that fail are tried again next time)
"""
import argparse
import gzip
import json
import re
import threading
import time

import numpy as np
import pandas as pd
import requests

from soar_meta import OUT_DIR, months

API = "https://www.lmsal.com/hek/hcr"
CACHE = OUT_DIR / "cache" / "iris"
FIRST_MONTH = pd.Timestamp("2013-07-01")   # IRIS launched 2013-06-27, first light 2013-07-17
LIMIT = 5000                         # a month has ~300-600 observations; reaching this means cut off
SJI_ABOUT = {1330: "C II 1330 (upper chromosphere / transition region)", 1400: "Si IV 1400 (transition region)",
             2796: "Mg II k 2796 (chromosphere)", 2832: "Mg II wing 2832 (upper photosphere)",
             1600: "1600 continuum (photosphere)", 5000: "5000 white light"}
MONTHS = {}                           # month -> DataFrame, in memory (also read by browse_web for the timeline)
_lock = threading.Lock()


def query_month(m, refresh=False, verbose=True):
    """All IRIS observations of month m as HCR returns them (list of dicts), cached raw."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"hcr_{m:%Y-%m}.json.gz"
    if path.is_file() and not refresh:
        with gzip.open(path, "rt") as f:
            return json.load(f)
    if verbose:
        print(f"  IRIS HCR: {m:%Y-%m} ...", end=" ", flush=True)
    end = m + pd.offsets.MonthBegin(1)
    for attempt in range(3):   # HCR is slow (~40 s a month) and sometimes times out
        try:
            r = requests.get(API, params=dict(cmd="search-events3", outputformat="json", instrument="IRIS",
                                              startTime=f"{m:%Y-%m-%dT%H:%M}", stopTime=f"{end:%Y-%m-%dT%H:%M}",
                                              limit=LIMIT), timeout=600)
            r.raise_for_status()
            ev = r.json().get("Events", [])
            break
        except (requests.RequestException, ValueError):
            if attempt == 2:
                raise
            time.sleep(30)
    if verbose:
        print(f"{len(ev)} observations")
    if len(ev) >= LIMIT:
        print(f"  WARNING: IRIS {m:%Y-%m}: {len(ev)} records, probably cut off; not cached")
        return ev
    tmp = path.with_suffix(".tmp")
    with gzip.open(tmp, "wt") as f:
        json.dump(ev, f)
    tmp.replace(path)
    return ev


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def _text(v):
    return str(v).strip().strip("'").strip() if v is not None else ""


def _row(e):
    """One observation: the useful fields, flat, plus the channel groups (list of dicts)."""
    groups = e.get("groups") or []
    ras = [g for g in groups if str(g.get("group_name")).strip("'") == "Raster"]
    sji = [g for g in groups if _num(g.get("wavelnth")) > 0 and "deconvolved" not in str(g.get("group_name"))]
    g = ras[0] if ras else {}
    channels = sorted({int(_num(s["wavelnth"])) for s in sji})
    sji_cad = [_num(s.get("cadence_avg_asrun")) for s in sji]
    return dict(
        event_key=e.get("eventKey"), obsid=_text(e.get("obsId")),
        start=pd.Timestamp(e["startTime"]), end=pd.Timestamp(e["stopTime"]),
        aim=_text(e.get("sciObjectives")), planners=_text(e.get("planners")), goal=_text(e.get("goal")),
        target=_text(e.get("target")), description=_text(e.get("descriptions")),
        xcen=_num(e.get("xCen")), ycen=_num(e.get("yCen")), roll=_num(e.get("roll_angle")),
        raster_w=_num(e.get("raster_fovx")), raster_h=_num(e.get("raster_fovy")),
        sji_w=_num(e.get("sji_fovx")), sji_h=_num(e.get("sji_fovy")),
        cover_w=_num(e.get("xFov")), cover_h=_num(e.get("yFov")),
        steps=_num(g.get("raster_numsteps")), step_size=_num(g.get("raster_stepsize")),
        step_cadence=_num(g.get("basic_cadence_avg_asrun")), raster_cadence=_num(g.get("cadence_avg_asrun")),
        n_rasters=_num(g.get("numrasters")), exptime=_num(g.get("max_exptime")),
        sji_channels=",".join(map(str, channels)),
        sji_cadence=float(np.nanmedian(sji_cad)) if sji_cad and np.isfinite(sji_cad).any() else np.nan,
        sji_images=int(np.nansum([_num(s.get("num_images")) for s in sji])) if sji else 0,
        lines=", ".join(dict.fromkeys(x.strip() for x in _text(e.get("spectral_window_names")).split(",") if x.strip())),
        size_mb=float(np.nansum([_num(x.get("comp_data_size_mb")) for x in groups])) if groups else np.nan,
        url=_text(e.get("parentUrl")), groups=groups)


def month_table(m, refresh=False, verbose=True):
    """The observations of month m as a DataFrame (in memory after the first call)."""
    m = pd.Timestamp(m).to_period("M").start_time
    with _lock:
        if m in MONTHS and not refresh:
            return MONTHS[m]
    ev = query_month(m, refresh, verbose)
    df = pd.DataFrame([_row(e) for e in ev])
    if len(df):
        df = df.sort_values("start").reset_index(drop=True)
    with _lock:
        MONTHS[m] = df
    return df


def observations(t0, t1, refresh=False, verbose=True):
    """Observations overlapping [t0, t1], one row each (see _row for the columns)."""
    t0, t1 = pd.Timestamp(t0), pd.Timestamp(t1)
    parts = [month_table(m, refresh, verbose) for m in months(t0 - pd.Timedelta("2D"), t1)]
    parts = [p for p in parts if len(p)]
    if not parts:
        return pd.DataFrame()
    df = pd.concat(parts, ignore_index=True).drop_duplicates("event_key")
    return df[(df.start <= t1) & (df.end >= t0)].sort_values("start").reset_index(drop=True)


def corners(r, which="raster"):
    """Corners (4 x 2, arcsec, Earth view) of the raster or SJI area of observation r,
    rotated by the roll (taken counter-clockwise on the sky; all 2025-10 rolls are 0)."""
    w, h = (r.raster_w, r.raster_h) if which == "raster" else (r.cover_w, r.cover_h)
    if not (np.isfinite(w) and np.isfinite(h)):
        return None
    a = np.radians(r.roll if np.isfinite(r.roll) else 0.0)
    rot = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    box = np.array([[-w, -h], [w, -h], [w, h], [-w, h]]) / 2
    return box @ rot.T + [r.xcen, r.ycen]


def is_sit_and_stare(r):
    return (np.isfinite(r.step_size) and abs(r.step_size) < 0.05) or "sit-and-stare" in r.goal.lower()


def raster_text(r):
    """The spectrograph programme in words."""
    if not np.isfinite(r.steps):
        return "no raster information"
    if is_sit_and_stare(r):   # HCR counts the exposures as "steps" (step size 0)
        out = [f"sit-and-stare (slit fixed, {r.raster_h:.0f}\" long), {int(r.steps)} exposures"]
        if np.isfinite(r.step_cadence):
            out.append(f"one every {r.step_cadence:.3g} s")
    else:
        out = [f"{int(r.steps)}-step raster, {abs(r.step_size):.2g}\" steps, {r.raster_w:.0f}\"x{r.raster_h:.0f}\""]
        if np.isfinite(r.step_cadence):
            out.append(f"{r.step_cadence:.3g} s per step")
        if np.isfinite(r.raster_cadence):
            out.append(f"one raster every {r.raster_cadence:.0f} s")
        if np.isfinite(r.n_rasters):
            out.append(f"{int(r.n_rasters)} raster{'s' if r.n_rasters != 1 else ''}")
    if np.isfinite(r.exptime):
        out.append(f"exposure {r.exptime:.2g} s")
    return ", ".join(out)


def sji_text(r):
    """The slit-jaw imaging in words."""
    if not r.sji_channels:
        return "no slit-jaw images"
    ch = [int(c) for c in r.sji_channels.split(",")]
    cad = f", each channel every {r.sji_cadence:.3g} s" if np.isfinite(r.sji_cadence) else ""
    return (f"SJI {' / '.join(map(str, ch))} Å{cad}, {r.sji_images} images, "
            f"{r.cover_w:.0f}\"x{r.cover_h:.0f}\" covered")


def fetch_all(refresh=False):
    """Cache every month from FIRST_MONTH to the current one (those already cached are skipped)."""
    ms = months(FIRST_MONTH, pd.Timestamp.now())
    todo = [m for m in ms if refresh or not (CACHE / f"hcr_{m:%Y-%m}.json.gz").is_file()]
    print(f"IRIS HCR: {len(ms)} months {ms[0]:%Y-%m} .. {ms[-1]:%Y-%m}, {len(ms) - len(todo)} in the cache, "
          f"{len(todo)} to fetch")
    failed, n = [], 0
    for k, m in enumerate(todo):
        print(f"[{k + 1}/{len(todo)}]", end=" ", flush=True)
        try:
            n += len(query_month(m, refresh=refresh))
        except Exception as e:
            print(f"  {m:%Y-%m} failed: {e}")
            failed.append(m)
    print(f"done: {n} observations fetched; " + (f"failed (run again): {', '.join(f'{m:%Y-%m}' for m in failed)}"
                                                 if failed else "all months cached"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("start", nargs="?")
    ap.add_argument("end", nargs="?")
    ap.add_argument("--all", action="store_true", help="cache every month since 2013-07 (not yet cached)")
    ap.add_argument("--pole", action="store_true", help='only |ycen| > 800" (polar pointings)')
    ap.add_argument("--long", action="store_true", help="also aim, planner, programme and channels")
    ap.add_argument("--refresh", action="store_true", help="query HCR again (months of the range)")
    args = ap.parse_args()
    if args.all:
        fetch_all(args.refresh)
        return
    if not (args.start and args.end):
        ap.error("give start and end (or --all)")
    d = observations(args.start, args.end, args.refresh)
    if args.pole and len(d):
        d = d[d.ycen.abs() > 800]
    if d.empty:
        print("No IRIS observations.")
        return
    print(f"{'start':16}  {'end':11}  {'OBSID':10}  {'xcen':>6} {'ycen':>6}  {'roll':>4}  {'raster':>9}  aim / planner")
    for _, r in d.iterrows():
        end = f"{r.end:%H:%M}" if r.end.date() == r.start.date() else f"{r.end:%m-%d %H:%M}"
        print(f"{r.start:%Y-%m-%d %H:%M}  {end:11}  {r.obsid:10}  {r.xcen:6.0f} {r.ycen:6.0f}  {r.roll:4.0f}  "
              f"{r.raster_w:4.0f}x{r.raster_h:<4.0f}  {r.aim[:45]} / {r.planners}")
        if args.long:
            print(f"{'':18}{r.goal}")
            print(f"{'':18}{raster_text(r)}")
            print(f"{'':18}{sji_text(r)}; level 2 {r.size_mb:.0f} MB; {r.url}")
    print(f"{len(d)} observations")


if __name__ == "__main__":
    main()
