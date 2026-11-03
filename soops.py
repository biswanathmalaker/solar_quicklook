#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
ESA's SOOP summary table (https://www.cosmos.esa.int/web/solar-orbiter/soops-summary,
saved by hand as soops-summary.html in this folder) as a table: one row per SOOP run, with
what it targeted, the coordinated observations (IRIS, Hinode, DKIST, others) and notes.

The HTML is parsed once into cache/soops/soops_summary.csv; it is parsed again only when
soops-summary.html is newer than the CSV (i.e. after you save a newer copy of the page).

    import soops
    s = soops.table()                                    # all runs
    soops.runs("R_SMALL_HRES_MCAD_Polar-Rotation", "2025-10-20")   # runs of a SOOP around a time

Run:  ./soops.py                       # list all runs
      ./soops.py Polar                 # runs whose SOOP name contains "Polar"
      ./soops.py --date 2025-10-20     # runs covering this date
"""
import argparse
import html
import json
import re

import pandas as pd
import requests

from soar_meta import OUT_DIR

HTML = OUT_DIR / "soops-summary.html"
CSV = OUT_DIR / "cache" / "soops" / "soops_summary.csv"
PAGES = OUT_DIR / "cache" / "soops" / "pages"     # one JSON per SOOP description page
PAGE_URL = "https://www.cosmos.esa.int/web/solar-orbiter/soop_pages/{}"
COLUMNS = ["soop", "description", "rsw", "ltp", "start", "end", "coordinators", "iris", "hinode", "dkist",
           "other_coordination", "status", "data_link", "quicklook", "notes"]


def _text(cell):
    t = html.unescape(re.sub(r"<[^>]+>", " ", cell))
    return re.sub(r"\s+", " ", t).strip()


def _links(cell):
    return " ".join(html.unescape(h) for h in re.findall(r'href="([^"]+)"', cell))


def parse(path=HTML):
    """soops-summary.html -> DataFrame, one row per SOOP run. soop_full is the full SOOP
    name as in the SOAR metadata (e.g. R_SMALL_HRES_MCAD_Polar-Rotation), taken from the
    observation-summary link; soop_page is the link to the SOOP's description page."""
    s = path.read_text(encoding="utf-8", errors="ignore")
    table = s[s.find("<table"):]
    rows = []
    for tr in re.findall(r"<tr.*?</tr>", table, flags=re.S)[1:]:
        cells = re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", tr, flags=re.S)
        if len(cells) < len(COLUMNS):
            continue
        row = {c: _text(v) for c, v in zip(COLUMNS, cells)}
        row["soop_page"] = _links(cells[0])
        row["data_link"] = _links(cells[12]) or row["data_link"]
        row["quicklook"] = _links(cells[13])
        row["other_coordination_links"] = _links(cells[10])
        m = re.search(r"soop=([A-Z]_[A-Z]+_[A-Z]+_[A-Z]+_[^&\"]+?)-(\d{8}T\d{4})", row["quicklook"])
        row["soop_full"] = m.group(1) if m else ""
        rows.append(row)
    df = pd.DataFrame(rows)
    for c in ("start", "end"):
        df[c] = pd.to_datetime(df[c].str.replace("*", "", regex=False), errors="coerce")
    return df.sort_values("start").reset_index(drop=True)


def page(soop, refresh=False):
    """What a SOOP is for, from its ESA description page (fetched once, kept in
    cache/soops/pages/): dict(soop, url, description (the "General description" text),
    instruments ([instrument, mode, comment] rows of the baseline observations table)).
    None if the page cannot be read (not cached, tried again next session)."""
    soop = str(soop).strip()
    if not soop or soop.lower() == "none":
        return None
    path = PAGES / f"{soop}.json"
    if path.is_file() and not refresh:
        return json.loads(path.read_text())
    url = PAGE_URL.format(soop.lower())
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
    except requests.RequestException:
        return None
    t = r.text
    i = t.find("General description")
    if i < 0:
        return None
    body = t[i:]
    end = re.search(r"Default SOOP duration|Observations requirement|Pointing requirements", body)
    desc = _text(re.sub(r"</p>|<br\s*/?>|</li>", " ", body[:end.start() if end else 4000]))
    desc = desc.replace("General description", "", 1).strip()
    instruments = []
    j = body.find("Observations requirement")
    if j >= 0:
        tab = re.search(r"<table.*?</table>", body[j:], flags=re.S)
        for tr in re.findall(r"<tr.*?</tr>", tab.group(0) if tab else "", flags=re.S):
            cells = [_text(c) for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", tr, flags=re.S)]
            if cells and cells[0] and cells[0].lower() != "instrument":
                instruments.append((cells + ["", "", ""])[:3])
    out = dict(soop=soop, url=url, description=desc, instruments=instruments)
    PAGES.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    return out


def table():
    """The SOOP summary as a DataFrame (from the CSV; re-parsed if the HTML is newer)."""
    if not CSV.is_file() or (HTML.is_file() and HTML.stat().st_mtime > CSV.stat().st_mtime):
        df = parse()
        CSV.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(CSV, index=False)
        print(f"  SOOP summary: {len(df)} runs parsed from {HTML.name} -> {CSV.relative_to(OUT_DIR)}")
        return df
    return pd.read_csv(CSV, parse_dates=["start", "end"], keep_default_na=False)


def runs(soop=None, t0=None, t1=None, margin="1D"):
    """Runs of a SOOP (full name as in SOAR, or part of the short name) overlapping
    [t0, t1] widened by margin (the table's times are planning times)."""
    df = table()
    if soop:
        full = df.soop_full.str.lower() == soop.lower()
        df = df[full] if full.any() else df[df.soop.str.contains(soop.split("_")[-1], case=False, regex=False)]
    if t0 is not None:
        t0 = pd.Timestamp(t0)
        t1 = pd.Timestamp(t1) if t1 is not None else t0
        w = pd.Timedelta(margin)
        df = df[(df.start <= t1 + w) & (df.end >= t0 - w)]
    return df


def describe(r, width=110):
    """A few lines about one run."""
    coord = [f"{k} {r[c]}" for k, c in (("IRIS:", "iris"), ("Hinode:", "hinode"), ("DKIST:", "dkist"),
                                         ("other:", "other_coordination")) if str(r[c]).strip() not in ("", "nan")]
    lines = [f"{r.soop}  {r.start:%Y-%m-%d %H:%M} .. {r.end:%Y-%m-%d %H:%M}  {r.rsw}  coord.: {r.coordinators}",
             f"    {r.description}" if str(r.description).strip() else "",
             f"    coordinated: {'; '.join(coord)}" if coord else "",
             f"    notes: {r.notes}" if str(r.notes).strip() not in ("", "nan") else ""]
    return "\n".join(line[:width] for line in lines if line)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("soop", nargs="?", help="part of the SOOP name (case insensitive)")
    ap.add_argument("--date", help="runs covering this date (±1 day)")
    args = ap.parse_args()
    df = runs(args.soop, args.date) if (args.soop or args.date) else table()
    for _, r in df.iterrows():
        print(describe(r, width=200))


if __name__ == "__main__":
    main()
