#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
What Solar Orbiter observed in a time range: the observing blocks (one pointing, no long
gap) of HRIEUV, HRI-Lya, PHI-HRT, FSI, PHI-FDT and SPICE, with where the high-resolution
telescopes pointed and whether Earth (SDO) sees the same region. Only metadata (SOAR TAP,
cached in cache/soar/), no FITS download.

Columns of the table
  block            block number (use it with plot_fov.py --block)
  start, dur       first file (UTC) and length of the block
  n, cad           number of files, median time between them
  soop             SOOP name without its R_/L_ prefix part (e.g. Polar-Rotation)
  SolO             Solar Orbiter: Stonyhurst lat, lon (= angle from Earth, + west), r (au)
  centre           Stonyhurst lat, lon of the FOV centre; "off" when the centre is off disk
  disk             fraction of the FOV on the disk (as seen from Solar Orbiter)
  lat range        latitudes covered by the on-disk part
  Earth            fraction of the on-disk part that Earth also sees, and mu seen from
                   Earth at the centre (0 = at Earth's limb)
  src              where the pointing comes from (fov / header / inferred, see soar_meta.py)

Run:  ./find_obs.py 2025-10-01 2025-10-31                       # hri174 and phi-hrt
      ./find_obs.py 2025-03-01 2025-04-30 --products hri174 phi-hrt fsi174
      ./find_obs.py 2025-01-01 2025-12-31 --soop Polar          # SOOP name contains "Polar"
      ./find_obs.py 2025-01-01 2025-12-31 --min-lat 60          # FOV reaches |lat| >= 60
      ./find_obs.py 2025-01-01 2025-12-31 --both hri174 phi-hrt # HRIEUV + PHI-HRT together
      ./find_obs.py ... --csv blocks.csv                        # also save the table
      ./find_obs.py ... --refresh                               # query SOAR again (new releases)

Below the tables: the runs of ESA's SOOP summary (soops.py) behind the blocks' SOOPs --
what was targeted, coordinated observations (IRIS, Hinode, DKIST, ...), notes.
"""
import argparse

import numpy as np
import pandas as pd

import soar_meta as sm
import soops


def short_soop(s):
    parts = str(s).split("_")
    return parts[-1] if len(parts) >= 5 else str(s)


def fmt_dur(td):
    s = td.total_seconds()
    return f"{s / 3600:.1f}h" if s >= 3600 else f"{s / 60:.0f}m"


def fmt_cad(s):
    if not np.isfinite(s):
        return "-"
    return f"{s:.0f}s" if s < 120 else f"{s / 60:.0f}m" if s < 7200 else f"{s / 3600:.1f}h"


def print_blocks(b):
    if b.empty:
        print("  (none)")
        return
    pointed = "lat" in b
    head = f"{'block':>5}  {'start':16}  {'dur':>5} {'n':>5} {'cad':>5}  {'soop':24}"
    if pointed:
        head += f"  {'SolO lat  lon    r':18}  {'centre lat  lon':15}  {'disk':>4}  {'lat range':11}  {'Earth  mu':9}  src"
    print(head)
    for _, r in b.iterrows():
        line = (f"{r.block:>5}  {r.start:%Y-%m-%d %H:%M}  {fmt_dur(r.end - r.start):>5} {r.n:>5} "
                f"{fmt_cad(r.cadence_s):>5}  {short_soop(r.soop)[:24]:24}")
        if pointed and np.isfinite(r.get("solo_lat", np.nan)):
            centre = f"{r.lat:+5.0f} {r.lon:+5.0f}" if np.isfinite(r.lat) else "  off"
            rng = f"{r.lat_min:+4.0f}..{r.lat_max:+3.0f}" if np.isfinite(r.lat_min) else "-"
            earth = f"{r.earth_vis:4.0%} {r.mu_earth:4.2f}" if np.isfinite(r.earth_vis) else "-"
            earth = earth.replace("nan", " -- ")
            line += (f"  {r.solo_lat:+5.1f} {r.solo_lon:+6.1f} {r.solo_r:4.2f}  {centre:15}  "
                     f"{r.on_disk:4.0%}  {rng:11}  {earth:9}  {r.src}")
        print(line)


def print_soop_runs(b):
    """The runs of ESA's SOOP summary (soops.py) behind the blocks' SOOPs: what was targeted,
    coordinated observations, notes."""
    shown = set()
    lines = []
    for soop, g in b[~b.soop.isin(["none", "not defined"])].groupby("soop"):
        for _, r in soops.runs(soop, g.start.min(), g.end.max(), margin="2h").iterrows():
            key = (r.soop_full or r.soop, r.start)
            if key not in shown:
                shown.add(key)
                lines.append((r.start, soops.describe(r, width=160)))
    if lines:
        print("\nSOOP runs (ESA SOOP summary, soops-summary.html):")
        for _, text in sorted(lines):
            print(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("start")
    ap.add_argument("end", help="dates are UTC; a bare date is its 00:00")
    ap.add_argument("--products", nargs="+", default=["hri174", "phi-hrt"], choices=list(sm.PRODUCTS))
    ap.add_argument("--soop", help="keep blocks whose SOOP name contains this (case insensitive)")
    ap.add_argument("--min-lat", type=float, help="keep pointed blocks whose on-disk FOV reaches |lat| >= this")
    ap.add_argument("--min-files", type=int, default=1, help="keep blocks with at least this many files")
    ap.add_argument("--both", nargs=2, metavar=("A", "B"),
                    help="only blocks of A overlapping (in time and on the Sun) a block of B, e.g. hri174 phi-hrt")
    ap.add_argument("--csv", help="save the blocks table to this file")
    ap.add_argument("--refresh", action="store_true", help="query SOAR again instead of the cache")
    ap.add_argument("--no-headers", action="store_true", help="skip the slow FITS-header queries")
    ap.add_argument("--no-soop-info", action="store_true", help="do not list the SOOP runs (ESA's SOOP summary)")
    args = ap.parse_args()

    products = list(dict.fromkeys(args.products + (args.both or [])))
    f = sm.files(args.start, args.end, products, refresh=args.refresh, headers=not args.no_headers)
    if f.empty:
        print("No files.")
        return
    b, f = sm.blocks(f)
    b = b[b.n >= args.min_files]
    if args.soop:
        b = b[b.soop.str.contains(args.soop, case=False, na=False)]
    if args.min_lat is not None and "lat_min" in b:
        reach = np.fmax(b.lat_min.abs(), b.lat_max.abs())
        b = b[~b["product"].isin(sm.POINTED) | (reach >= args.min_lat)]
    pairs = None
    if args.both:
        a_, b_ = args.both
        pairs = sm.overlap(b[b["product"] == a_], b[b["product"] == b_])
        keep = set(pairs.block_a) | set(pairs.block_b) if len(pairs) else set()
        b = b[b.block.isin(keep)]

    for product in products:
        sel = b[b["product"] == product].drop(columns=[c for c in b if c.startswith("c") and c[1:2].isdigit()])
        print(f"\n{product}  ({', '.join(sm.PRODUCTS[product]['desc'])}): {len(sel)} blocks, "
              f"{int(sel.n.sum()) if len(sel) else 0} files")
        print_blocks(sel.dropna(axis=1, how="all"))
    if pairs is not None and len(pairs):
        print(f"\n{args.both[0]} block / {args.both[1]} block, common time, fraction of the "
              f"{args.both[0]} FOV inside the {args.both[1]} FOV:")
        for _, p in pairs.iterrows():
            print(f"  {p.block_a:>5} / {p.block_b:<5}  {p.start:%Y-%m-%d %H:%M} .. {p.end:%Y-%m-%d %H:%M}  "
                  f"{p.frac:4.0%}")
    if not args.no_soop_info:
        print_soop_runs(b)
    if args.csv:
        b.to_csv(args.csv, index=False)
        print(f"\nSaved {len(b)} blocks -> {args.csv}")


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    main()
