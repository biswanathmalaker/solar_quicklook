#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
DKIST datasets (VBI, ViSP, DL-NIRSP, Cryo-NIRSP, VTF) in a time range, from the DKIST
Data Center search API (https://api.dkistdc.nso.edu/datasets/v1, no login needed;
embargoed datasets are listed too, flagged). Each dataset has a bounding box in
helioprojective arcsec as seen from Earth, at its start. Queried by month and kept in
cache/dkist/ like the SOAR metadata (soar_meta.py): once, again only with --refresh.

    import dkist_meta as dm
    d = dm.datasets("2025-10-01", "2025-10-31")     # one row per dataset

Run:  ./dkist_meta.py 2025-10-01 2025-10-31          # table of the observations
      ./dkist_meta.py --all                          # cache every month since 2020-01 not yet cached
"""
import argparse
import re

import numpy as np
import pandas as pd
import requests

from soar_meta import OUT_DIR, months

API = "https://api.dkistdc.nso.edu/datasets/v1"
CACHE = OUT_DIR / "cache" / "dkist"
PAGE = 500
FIRST_MONTH = pd.Timestamp("2020-01-01")   # DKIST first light 2019-12; science data from 2021-22
KEEP = ["datasetId", "instrumentName", "startTime", "endTime", "boundingBox", "primaryProposalId",
        "isEmbargoed", "embargoEndDate", "spectralLines", "wavelengthMin", "wavelengthMax", "frameCount",
        "averageDatasetSpatialSampling", "qualityAverageFriedParameter", "browseMovieUrl",
        "primaryExperimentId", "experimentDescription", "targetTypes", "datasetSize"]   # text: search_meta.py


def query_month(m, refresh=False, verbose=True):
    """All datasets starting in month m (pages of PAGE), cached in cache/dkist/."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"datasets_{m:%Y-%m}.csv.gz"
    if path.is_file() and not refresh:   # queried once; again only with refresh
        return pd.read_csv(path)
    if verbose:
        print(f"  DKIST: {m:%Y-%m} ...", end=" ", flush=True)
    end = m + pd.offsets.MonthBegin(1)
    rows, page, total = [], 1, None
    while total is None or len(rows) < total:
        r = requests.get(API, params=dict(startTimeMin=f"{m:%Y-%m-%dT%H:%M:%S}", startTimeMax=f"{end:%Y-%m-%dT%H:%M:%S}",
                                          pageSize=PAGE, pageNumber=page), timeout=120)
        r.raise_for_status()
        d = r.json()
        total = d["recordCount"]
        if not d["searchResults"]:
            break
        rows += d["searchResults"]
        page += 1
    df = pd.DataFrame(rows, columns=None if rows else KEEP)
    df = df[[c for c in KEEP if c in df]]
    for c in ("spectralLines", "targetTypes"):
        if c in df:
            df[c] = df[c].apply(lambda v: "; ".join(v) if isinstance(v, list) else v)
    if verbose:
        print(f"{len(df)} datasets")
    if total is not None and len(df) < total:
        print(f"  WARNING: DKIST {m:%Y-%m}: {len(df)} of {total} datasets; not cached")
        return df
    tmp = path.with_suffix(".tmp")
    df.to_csv(tmp, index=False, compression="gzip")
    tmp.replace(path)
    return df


def parse_box(text):
    """'(x1,y1),(x2,y2)' (arcsec, two opposite corners) -> x0, y0, x1, y1 (min, max)."""
    v = [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", str(text))]
    if len(v) != 4:
        return (np.nan,) * 4
    return min(v[0], v[2]), min(v[1], v[3]), max(v[0], v[2]), max(v[1], v[3])


def datasets(t0, t1, refresh=False):
    """Datasets overlapping [t0, t1]: datasetId, instrument, start, end, x0, y0, x1, y1
    (bounding box, arcsec, Earth view), proposal, embargoed, lines, ..."""
    t0, t1 = pd.Timestamp(t0), pd.Timestamp(t1)
    # a dataset can start in the month before t0 and still overlap it
    ms = months(t0 - pd.Timedelta("1D"), t1)
    df = pd.concat([query_month(m, refresh) for m in ms], ignore_index=True)
    if df.empty:
        return df
    df["start"] = pd.to_datetime(df.startTime)
    df["end"] = pd.to_datetime(df.endTime)
    df = df[(df.start <= t1) & (df.end >= t0)].copy()
    df[["x0", "y0", "x1", "y1"]] = [parse_box(b) for b in df.boundingBox] if len(df) else np.empty((0, 4))
    df = df.rename(columns={"instrumentName": "instrument", "primaryProposalId": "proposal",
                            "isEmbargoed": "embargoed", "spectralLines": "lines"})
    return df.drop_duplicates("datasetId").sort_values("start").reset_index(drop=True)


def fetch_all(refresh=False):
    """Cache every month from FIRST_MONTH to now (cached months are skipped unless refresh)."""
    ms = months(FIRST_MONTH, pd.Timestamp.now())
    todo = [m for m in ms if refresh or not (CACHE / f"datasets_{m:%Y-%m}.csv.gz").is_file()]
    print(f"DKIST: {len(ms)} months {ms[0]:%Y-%m} .. {ms[-1]:%Y-%m}, {len(ms) - len(todo)} in the cache, {len(todo)} to fetch")
    failed, n = [], 0
    for m in todo:
        try:
            n += len(query_month(m, refresh=refresh))
        except Exception as e:
            print(f"  {m:%Y-%m} failed: {e}")
            failed.append(m)
    print(f"done: {n} datasets fetched; " + (f"failed (run again): {', '.join(f'{m:%Y-%m}' for m in failed)}"
                                             if failed else "all months cached"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("start", nargs="?")
    ap.add_argument("end", nargs="?")
    ap.add_argument("--all", action="store_true", help="cache every month since 2020-01 (not yet cached)")
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()
    if args.all:
        fetch_all(args.refresh)
        return
    if not (args.start and args.end):
        ap.error("give start and end (or --all)")
    d = datasets(args.start, args.end, args.refresh)
    if d.empty:
        print("No DKIST datasets.")
        return
    print(f"{'start':16}  {'end':5}  {'instrument':10}  {'x0':>7} {'y0':>7} {'x1':>7} {'y1':>7}  {'proposal':10}  lines")
    for _, r in d.iterrows():
        print(f"{r.start:%Y-%m-%d %H:%M}  {r.end:%H:%M}  {r.instrument:10}  {r.x0:7.0f} {r.y0:7.0f} {r.x1:7.0f} {r.y1:7.0f}"
              f"  {r.proposal:10}  {str(r.lines)[:40]}{'  (embargoed)' if r.embargoed else ''}")


if __name__ == "__main__":
    main()
