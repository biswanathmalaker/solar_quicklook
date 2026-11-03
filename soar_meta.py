#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
Solar Orbiter archive (SOAR) metadata, without downloading any FITS file: which EUI / PHI /
SPICE files exist in a time range, where each high-resolution image pointed, and how the
files group into observing blocks (one pointing, no long gap).

Where the information comes from (SOAR TAP service, queried by month and kept in
cache/soar/: queried once, again only with refresh=True / --refresh):
  v_sc_data_item   every L2 file: name, start / end, SOOP. The complete list.
  v_eui_hri_fov    HRI (EUV and Lya) and PHI-HRT pointing: two opposite corners of the
  v_phi_hrt_fov    field of view (pixels (0,0) and (N-1,N-1), helioprojective as seen
                   from Solar Orbiter, arcsec), the roll angle CROTA and the positions of
                   Solar Orbiter, Earth, ... Fast, but not complete (e.g. 56% of the
                   HRIEUV files of 2025-10 are missing).
  v_eui_sc_fits    the FITS headers (CRVAL, CRPIX, CDELT, PC, NAXIS, observer): the
  v_phi_sc_fits    pointing of the files missing from the FOV tables. Slow (~1 min a
                   month for EUI) and also not complete for PHI (no 2024-10, 2025-09/10),
                   so queried only for months with files still missing.
Files with neither take the pointing of the nearest file of the same block ("inferred").

Use from Python:
    import soar_meta as sm
    files = sm.files("2025-10-01", "2025-10-31", ["hri174", "phi-hrt"])
    blocks = sm.blocks(files)
or from the command line: ./find_obs.py (see there).
"""
import io
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

OUT_DIR = Path(__file__).resolve().parent
CACHE = OUT_DIR / "cache" / "soar"
TAP = "https://soar.esac.esa.int/soar-sl-tap/tap/sync"
DATA_URL = "https://soar.esac.esa.int/soar-sl-tap/data"
EUI_URL = "https://data.observatory.be/science/solo/eui"   # ROB mirror of the EUI data
AU_KM = 1.495978707e8
RSUN_KM = 6.957e5

# name -> SOAR descriptors (L2), FOV table, header table and the file-name prefix of the
# product in the header table. Neither: full disk (FSI, FDT). SPICE: header table only (the
# raster: NAXIS1 slit positions x CDELT1 by NAXIS2 pixels along the slit x CDELT2).
PRODUCTS = {
    "hri174": dict(desc=("eui-hrieuv174-image", "eui-hrieuv174-image-short"),
                   fov="v_eui_hri_fov", hdr="v_eui_sc_fits", prefix="solo_L2_eui-hrieuv174"),
    # HRIEUV with the filter wheel open: still the 174 A passband (set by the mirrors), more
    # light; used e.g. for 10-s cadence runs (2025-10-20 12:30-14:30, polar coronal hole)
    "hriopn": dict(desc=("eui-hrieuvopn-image", "eui-hrieuvopn-image-short"),
                   fov="v_eui_hri_fov", hdr="v_eui_sc_fits", prefix="solo_L2_eui-hrieuvopn"),
    "hrilya": dict(desc=("eui-hrilya1216-image",),
                   fov="v_eui_hri_fov", hdr="v_eui_sc_fits", prefix="solo_L2_eui-hrilya1216"),
    "phi-hrt": dict(desc=("phi-hrt-blos",),   # one file per inversion; bazi, binc, ... share it
                    fov="v_phi_hrt_fov", hdr="v_phi_sc_fits", prefix="solo_L2_phi-hrt-blos"),
    "fsi174": dict(desc=("eui-fsi174-image",)),
    "fsi304": dict(desc=("eui-fsi304-image",)),
    "phi-fdt": dict(desc=("phi-fdt-blos",)),
    "spice": dict(desc=("spice-n-ras", "spice-n-sit", "spice-n-exp", "spice-w-ras", "spice-w-sit", "spice-w-exp"),
                  hdr="v_spi_sc_fits", prefix="solo_L2_spice"),   # no FOV table: headers only
}
POINTED = [p for p, v in PRODUCTS.items() if "hdr" in v]
# a new block starts after a gap longer than this, a change of SOOP, or a re-pointing
# (long enough to keep e.g. the hourly PHI-HRT maps of a Polar-Rotation run together)
BLOCK_GAP = pd.Timedelta("3h")
REPOINT_ARCSEC = 60
CHUNK_ROWS = 4000    # months with more rows are queried day by day


# ---------------------------------------------------------------- TAP and the cache

def tap(adql, timeout=900, retries=2):
    """Run an ADQL query on SOAR; returns a DataFrame (CSV)."""
    for attempt in range(retries + 1):
        try:
            r = requests.get(TAP, params=dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv",
                                              MAXREC=10_000_000, QUERY=adql), timeout=timeout)
            r.raise_for_status()
            break
        except requests.RequestException:
            if attempt == retries:
                raise
            time.sleep(5)
    if r.text.lstrip().startswith("<"):   # errors come back as a VOTable
        msg = re.search(r'value="ERROR">(.*?)</INFO>', r.text, re.S)
        raise RuntimeError(f"SOAR TAP: {msg.group(1).strip() if msg else r.text[:500]}")
    return pd.read_csv(io.StringIO(r.text), low_memory=False)


def months(t0, t1):
    """Month starts covering [t0, t1]."""
    return pd.date_range(pd.Timestamp(t0).to_period("M").start_time, pd.Timestamp(t1), freq="MS")


def fetch_verified(table, cols, where, tries=4):
    """SELECT cols FROM soar.table WHERE where, checked against SELECT COUNT(*) with the same
    WHERE: long responses are sometimes cut off without an error (rows come in no
    particular order, so the loss is scattered). Returns (rows, expected count)."""
    n = int(tap(f"SELECT COUNT(*) AS n FROM soar.{table} WHERE {where}").n.iloc[0])
    df = pd.DataFrame()
    for _ in range(tries):
        if n == 0:
            break
        df = tap(f"SELECT {cols} FROM soar.{table} WHERE {where}")
        if len(df) >= n:
            break
    return df, n


def cached(name, month, table, cols, where, time_col, refresh=False, verbose=True):
    """SELECT cols FROM soar.table WHERE where AND time_col in the month, from
    cache/soar/{name}_{YYYY-MM}.csv.gz if it is there: a month is queried once and then
    only again with refresh=True (--refresh), e.g. to get data released since (the file
    date says when it was queried). Months with more than CHUNK_ROWS rows are
    queried day by day (big responses are cut off far more often), every piece checked
    against its COUNT(*) (see fetch_verified). An incomplete month is not cached."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{name}_{month:%Y-%m}.csv.gz"
    if path.is_file() and not refresh:
        try:
            return pd.read_csv(path, low_memory=False)
        except pd.errors.EmptyDataError:   # a month with no rows
            return pd.DataFrame()
    if verbose:
        print(f"  SOAR: {name} {month:%Y-%m} ...", end=" ", flush=True)
    t = time.time()
    end = month + pd.offsets.MonthBegin(1)
    span = lambda a, b: f"{where} AND {time_col} >= '{a:%Y-%m-%dT%H:%M:%S}' AND {time_col} < '{b:%Y-%m-%dT%H:%M:%S}'"
    n = int(tap(f"SELECT COUNT(*) AS n FROM soar.{table} WHERE {span(month, end)}").n.iloc[0])
    if n <= CHUNK_ROWS:
        df, n = fetch_verified(table, cols, span(month, end))
        got, pieces = len(df), [df]
    else:
        pieces, got = [], 0
        for day in pd.date_range(month, end, freq="D", inclusive="left"):
            d, nd = fetch_verified(table, cols, span(day, day + pd.Timedelta("1D")))
            pieces.append(d)
            got += min(len(d), nd)
            if len(d) < nd and verbose:
                print(f"[{day:%d}: {len(d)}/{nd}]", end=" ", flush=True)
        n = max(n, got)
    df = pd.concat([p for p in pieces if len(p)], ignore_index=True) if any(len(p) for p in pieces) else pieces[0]
    if verbose:
        print(f"{len(df)} rows, {time.time() - t:.0f} s")
    if len(df) < n:
        print(f"  WARNING: {name} {month:%Y-%m}: {len(df)} of {n} rows; not cached, run again")
        return df
    tmp = path.with_suffix(".tmp")
    df.to_csv(tmp, index=False, compression="gzip")
    tmp.replace(path)
    return df



def item_id(names):
    """File name -> data item id (the name without _Vnn... and .fits)."""
    return names.str.replace(r"_V\d+.*$", "", regex=True)


# ---------------------------------------------------------------- the three sources

def query_items(product, m, refresh=False):
    desc = ", ".join(f"'{d}'" for d in PRODUCTS[product]["desc"])
    df = cached(f"items_{product}", m, "v_sc_data_item",
                "data_item_id, filename, descriptor, begin_time, end_time, soop_name, filesize",
                f"level='L2' AND descriptor IN ({desc})", "begin_time", refresh)
    if "filename" not in df:   # a month without files of this product
        return pd.DataFrame(columns=["data_item_id", "filename", "descriptor", "begin_time", "end_time",
                                     "soop_name", "filesize"])
    # several versions of a file: keep the last
    return df.sort_values("filename").drop_duplicates("data_item_id", keep="last")


FOV_COLS = ("filename, start_time, crota, fov_solo_bot_left_arcsec_tx, fov_solo_bot_left_arcsec_ty, "
            "fov_solo_top_right_arcsec_tx, fov_solo_top_right_arcsec_ty, "
            "solo_pos_lat_deg, solo_pos_lon_deg, solo_pos_radius_au, earth_pos_lat_deg, earth_pos_radius_au")


def query_fov(product, m, refresh=False):
    p = PRODUCTS[product]
    desc = "eui-hrieuv174-image" if product == "hri174" else p["desc"][0]   # short exposures included
    df = cached(f"fov_{product}", m, p["fov"], FOV_COLS,
                f"proc_level='L2' AND descriptor='{desc}'", "start_time", refresh)
    if "filename" not in df:   # no rows (e.g. hriopn: not in the FOV table)
        return pd.DataFrame(columns=["filename", "data_item_id", "crota"])
    df["data_item_id"] = item_id(df.filename)
    return df.sort_values("filename").drop_duplicates("data_item_id", keep="last")


def query_headers(product, m, refresh=False):
    p = PRODUCTS[product]
    df = cached(f"hdr_{product}", m, p["hdr"],
                "filename, dimension_index, num_pixels_per_axis, cr_val1, cr_val2, cr_pix1, cr_pix2, "
                "c_delt1, c_delt2, pc1_1, pc1_2, pc2_1, pc2_2, distance_sun_obsv, hg_latitude_obsv, hg_longitud_obsv",
                f"filename LIKE '{p['prefix']}%'", "date_average", refresh)
    if df.empty:
        return df
    # one row per image axis: NAXIS1 / NAXIS2 from dimension_index 1 / 2
    naxis = df.pivot_table(index="filename", columns="dimension_index", values="num_pixels_per_axis", aggfunc="first")
    df = df.drop_duplicates("filename").set_index("filename")
    df["naxis1"], df["naxis2"] = naxis.get(1), naxis.get(2, naxis.get(1))
    df = df.reset_index()
    df["data_item_id"] = item_id(df.filename)
    return df.sort_values("filename").drop_duplicates("data_item_id", keep="last")


def spice_modes(m, refresh=False):
    """SPICE files of month m: data_item_id, study (OBSERVATION_MODE, e.g.
    SCI_POLAR_SS_SL04_59.6S_FF), exposure (s), telapse (s), nx (slit positions; 1 for
    sit-and-stare). From the header table; the spectral windows of a study are in
    spice_meta.py."""
    df = cached("spicemode", m, "v_spi_sc_fits",
                "filename, observation_mode, exposure, elapse_time, num_pixels_per_axis",
                "dimension_index=1 AND filename LIKE 'solo_L2_spice%'", "date_average", refresh)
    if df.empty:
        return df
    df = df.rename(columns={"observation_mode": "study", "elapse_time": "telapse", "num_pixels_per_axis": "nx"})
    df["data_item_id"] = item_id(df.filename)
    return df.sort_values("filename").drop_duplicates("data_item_id", keep="last").drop(columns="filename")


# ---------------------------------------------------------------- corners

def corners_from_fov(df):
    """FOV table -> the four corners (N, 4, 2) in arcsec: the table gives pixels (0,0) and
    (N-1,N-1); the other two follow from the roll angle (checked against the WCS of a
    PHI-HRT file: same corners to 0.1")."""
    p0 = df[["fov_solo_bot_left_arcsec_tx", "fov_solo_bot_left_arcsec_ty"]].to_numpy(float)
    p2 = df[["fov_solo_top_right_arcsec_tx", "fov_solo_top_right_arcsec_ty"]].to_numpy(float)
    t = np.radians(df.crota.to_numpy(float))
    c, s = np.cos(t), np.sin(t)
    d = p2 - p0
    # diagonal in image axes: rotate back by CROTA
    d_x = c * d[:, 0] + s * d[:, 1]
    d_y = -s * d[:, 0] + c * d[:, 1]
    p1 = p0 + np.c_[c * d_x, s * d_x]
    p3 = p0 + np.c_[-s * d_y, c * d_y]
    return np.stack([p0, p1, p2, p3], axis=1)


def corners_from_header(df):
    """Header table -> the four corners (N, 4, 2) in arcsec (pixel centres of the corner
    pixels, linear WCS: fine over a few thousand arcsec)."""
    n1, n2 = df.naxis1.to_numpy(float), df.naxis2.to_numpy(float)
    out = []
    for p1, p2 in ((1, 1), (n1, 1), (n1, n2), (1, n2)):   # FITS pixels are 1-based
        dx = (p1 - df.cr_pix1.to_numpy(float))
        dy = (p2 - df.cr_pix2.to_numpy(float))
        tx = df.cr_val1 + df.c_delt1 * (df.pc1_1 * dx + df.pc1_2 * dy)
        ty = df.cr_val2 + df.c_delt2 * (df.pc2_1 * dx + df.pc2_2 * dy)
        out.append(np.c_[tx, ty])
    return np.stack(out, axis=1)


CORNER_COLS = [f"c{k}{a}" for k in range(4) for a in "xy"]


# ---------------------------------------------------------------- files

def files(t0, t1, products=("hri174", "phi-hrt"), refresh=False, headers=True):
    """All L2 files of the products in [t0, t1], one row per file, time ordered:
    product, data_item_id, filename, descriptor, start, end, soop, filesize and, for the
    pointed products (hri174, hrilya, phi-hrt, spice), the FOV corners c0x..c3y (arcsec, seen
    from Solar Orbiter), Solar Orbiter's position solo_lat, solo_lon (Stonyhurst, deg),
    solo_r (au), and pointing_src (fov / header / inferred / none).
    headers=False skips the slow header queries (more "inferred" pointings)."""
    t0, t1 = pd.Timestamp(t0), pd.Timestamp(t1)
    out = []
    for product in products:
        items = pd.concat([query_items(product, m, refresh) for m in months(t0, t1)], ignore_index=True)
        items["start"] = pd.to_datetime(items.begin_time)
        items["end"] = pd.to_datetime(items.end_time)
        items = items[(items.start >= t0) & (items.start <= t1)].copy()
        items["product"] = product
        items["soop"] = items.soop_name.fillna("none")
        if product in POINTED and len(items):
            items = add_pointing(items, product, t0, t1, refresh, headers)
        out.append(items)
    df = pd.concat(out, ignore_index=True) if out else pd.DataFrame()
    keep = ["product", "data_item_id", "filename", "descriptor", "start", "end", "soop", "filesize"]
    extra = [c for c in CORNER_COLS + ["solo_lat", "solo_lon", "solo_r", "pointing_src"] if c in df]
    return df[keep + extra].sort_values(["product", "start"]).reset_index(drop=True)


def add_pointing(items, product, t0, t1, refresh, headers):
    items = items.set_index("data_item_id")
    for c in CORNER_COLS + ["solo_lat", "solo_lon", "solo_r"]:
        items[c] = np.nan
    items["pointing_src"] = "none"
    # 1. FOV table
    if PRODUCTS[product].get("fov"):
        fov = pd.concat([query_fov(product, m, refresh) for m in months(t0, t1)], ignore_index=True)
        fov = fov[fov.data_item_id.isin(items.index)].set_index("data_item_id")
    else:
        fov = pd.DataFrame()
    if len(fov):
        items.loc[fov.index, CORNER_COLS] = corners_from_fov(fov).reshape(len(fov), 8)
        items.loc[fov.index, ["solo_lat", "solo_lon", "solo_r"]] = \
            fov[["solo_pos_lat_deg", "solo_pos_lon_deg", "solo_pos_radius_au"]].to_numpy()
        items.loc[fov.index, "pointing_src"] = "fov"
    # 2. headers, only for the months with files still missing
    if headers:
        missing = items[items.pointing_src == "none"]
        for m in months(t0, t1):
            if not ((missing.start >= m) & (missing.start < m + pd.offsets.MonthBegin(1))).any():
                continue
            hdr = query_headers(product, m, refresh)
            if hdr.empty:
                continue
            hdr = hdr[hdr.data_item_id.isin(missing.index)].set_index("data_item_id")
            if len(hdr):
                items.loc[hdr.index, CORNER_COLS] = corners_from_header(hdr).reshape(len(hdr), 8)
                items.loc[hdr.index, "solo_lat"] = hdr.hg_latitude_obsv
                items.loc[hdr.index, "solo_lon"] = hdr.hg_longitud_obsv
                items.loc[hdr.index, "solo_r"] = hdr.distance_sun_obsv / (AU_KM * 1e3)
                items.loc[hdr.index, "pointing_src"] = "header"
    return items.reset_index()


# ---------------------------------------------------------------- geometry

def hpc_to_hgs(tx, ty, lat_obs, lon_obs, r_au):
    """Helioprojective (arcsec, arrays) seen from an observer at Stonyhurst lat_obs, lon_obs
    (deg), r_au -> Stonyhurst lon, lat (deg) on the surface; NaN off disk. Closed form
    (Thompson 2006, A&A 449, 791), as in plot_ref_file.py."""
    D, R = r_au * AU_KM, RSUN_KM
    B0, L0 = np.radians(lat_obs), np.radians(lon_obs)
    tx, ty = np.radians(np.asarray(tx) / 3600), np.radians(np.asarray(ty) / 3600)
    q = D * np.cos(ty) * np.cos(tx)
    disc = q ** 2 - D ** 2 + R ** 2
    with np.errstate(invalid="ignore"):
        d = q - np.sqrt(disc)
        x, y, z = d * np.cos(ty) * np.sin(tx), d * np.sin(ty), D - d * np.cos(ty) * np.cos(tx)
        lat = np.arcsin(np.clip((y * np.cos(B0) + z * np.sin(B0)) / R, -1, 1))
        lon = L0 + np.arctan2(x, z * np.cos(B0) - y * np.sin(B0))
    lon, lat = (np.degrees(lon) + 180) % 360 - 180, np.degrees(lat)
    off = disc < 0
    return np.where(off, np.nan, lon), np.where(off, np.nan, lat)


def cos_angle(lon, lat, lon_obs, lat_obs):
    """Cosine of the angle at the Sun's centre between surface points and an observer."""
    lon, lat, lon_obs, lat_obs = map(np.radians, (lon, lat, lon_obs, lat_obs))
    return np.sin(lat) * np.sin(lat_obs) + np.cos(lat) * np.cos(lat_obs) * np.cos(lon - lon_obs)


def fov_grid(corners, n=15):
    """n x n points spread over a FOV (4, 2) -> (n*n, 2)."""
    u = (np.arange(n) + 0.5) / n
    a, b = np.meshgrid(u, u)
    p0, p1, p3 = corners[0], corners[1], corners[3]
    return p0 + a.ravel()[:, None] * (p1 - p0) + b.ravel()[:, None] * (p3 - p0)


def earth_position(t):
    """Stonyhurst lat (deg) and distance (au) of Earth at times t (lon is 0 by definition)."""
    from sunpy.coordinates import get_earth
    e = get_earth(pd.DatetimeIndex(t).to_pydatetime())
    return np.atleast_1d(e.lat.to_value("deg")), np.atleast_1d(e.radius.to_value("au"))


def describe_fov(corners, lat_obs, lon_obs, r_au, earth_lat):
    """Where a FOV looks, from its corners (4, 2) and Solar Orbiter's position: centre
    tx, ty, the Stonyhurst lon, lat of the centre (NaN off disk), the fraction of the FOV on
    the disk, the latitude range of its on-disk part, and the fraction of that part which
    Earth (SDO) also sees, plus mu seen from Earth at the centre."""
    centre = corners.mean(axis=0)
    lon_c, lat_c = hpc_to_hgs(centre[0], centre[1], lat_obs, lon_obs, r_au)
    g = fov_grid(corners)
    lon, lat = hpc_to_hgs(g[:, 0], g[:, 1], lat_obs, lon_obs, r_au)
    on = np.isfinite(lat)
    out = dict(tx=centre[0], ty=centre[1], lon=float(lon_c), lat=float(lat_c), on_disk=on.mean(),
               lat_min=np.nan, lat_max=np.nan, earth_vis=np.nan, mu_earth=np.nan)
    if on.any():
        out["lat_min"], out["lat_max"] = np.nanmin(lat), np.nanmax(lat)
        ce = cos_angle(lon[on], lat[on], 0.0, earth_lat)
        out["earth_vis"] = (ce > 0).mean()
        if np.isfinite(lat_c):
            out["mu_earth"] = float(cos_angle(lon_c, lat_c, 0.0, earth_lat))
    return out


# ---------------------------------------------------------------- blocks

def blocks(df, gap=BLOCK_GAP, repoint=REPOINT_ARCSEC):
    """Group the files of each product into observing blocks; files without pointing take
    that of the nearest file in their block ("inferred"). Returns (blocks, files with a
    block column). One row per block: product, block, start, end, n, cadence_s (median),
    soop, descriptors, solo_lat, solo_lon, solo_r, the corners, pointing_src counts, and
    the describe_fov() columns."""
    df = df.sort_values(["product", "start"]).reset_index(drop=True)
    df["block"] = -1
    rows, next_id = [], 0
    for product, f in df.groupby("product", sort=False):
        idx = f.index.to_numpy()
        t = f.start.to_numpy()
        new = np.r_[True, (t[1:] - f.end.to_numpy()[:-1]) > gap.to_timedelta64()]
        new |= np.r_[True, f.soop.to_numpy()[1:] != f.soop.to_numpy()[:-1]]
        if product in POINTED:
            # re-pointing: the centre of a file with pointing moved from that of the previous
            # file with pointing in the same segment (no gap, same SOOP) by more than repoint
            # and half the FOV's longer side (small offsets, e.g. the two alternating PHI-HRT pointings of a
            # Nanoflares run, stay in one block; mosaic tiles do not). Files without pointing
            # neither start a block nor are compared.
            seg = np.cumsum(new)
            cx = f[["c0x", "c1x", "c2x", "c3x"]].mean(axis=1).to_numpy()
            cy = f[["c0y", "c1y", "c2y", "c3y"]].mean(axis=1).to_numpy()
            side = np.fmax(np.hypot(f.c1x - f.c0x, f.c1y - f.c0y), np.hypot(f.c3x - f.c0x, f.c3y - f.c0y)).to_numpy()
            k = np.flatnonzero(np.isfinite(cx))   # files with pointing
            if len(k) > 1:
                jump = np.hypot(np.diff(cx[k]), np.diff(cy[k])) > np.fmax(repoint, 0.5 * side[k[1:]])
                jump &= seg[k[1:]] == seg[k[:-1]]
                new[k[1:][jump]] = True
        bid = np.cumsum(new) - 1 + next_id
        next_id = bid[-1] + 1 if len(bid) else next_id
        df.loc[idx, "block"] = bid
    if "c0x" in df:   # fill missing pointing within each block from the nearest file
        has = df.pointing_src.isin(["fov", "header"])
        cols = CORNER_COLS + ["solo_lat", "solo_lon", "solo_r"]
        for b, f in df[df["product"].isin(POINTED)].groupby("block"):
            ok = has[f.index]
            if ok.all() or not ok.any():
                continue
            src = f[ok]
            for i in f.index[~ok]:
                j = src.index[np.argmin(np.abs((src.start - df.at[i, "start"]).dt.total_seconds()))]
                df.loc[i, cols] = df.loc[j, cols].to_numpy()
                df.at[i, "pointing_src"] = "inferred"
    groups = list(df.groupby("block"))
    mids = [f.start.iloc[len(f) // 2] for _, f in groups]
    earth_lat = earth_position(mids)[0] if mids else []
    for (b, f), e_lat in zip(groups, earth_lat):
        product = f["product"].iloc[0]
        dt = f.start.diff().dt.total_seconds().dropna()
        row = dict(product=product, block=b, start=f.start.min(), end=f.end.max(), n=len(f),
                   cadence_s=float(dt.median()) if len(dt) else np.nan, soop=f.soop.iloc[0],
                   descriptors=",".join(sorted(f.descriptor.unique())))
        if product in POINTED and f[CORNER_COLS].notna().all(axis=1).any():
            g = f.dropna(subset=CORNER_COLS)
            corners = np.median(g[CORNER_COLS].to_numpy(float), axis=0).reshape(4, 2)
            lat_o, lon_o, r_o = (float(g[c].median()) for c in ("solo_lat", "solo_lon", "solo_r"))
            row.update(solo_lat=lat_o, solo_lon=lon_o, solo_r=r_o,
                       **dict(zip(CORNER_COLS, corners.ravel())),
                       src=",".join(f"{k}:{v}" for k, v in f.pointing_src.value_counts().items()),
                       **describe_fov(corners, lat_o, lon_o, r_o, e_lat))
        rows.append(row)
    return pd.DataFrame(rows), df


def overlap(blocks_a, blocks_b, min_frac=0.0):
    """Pairs of blocks (pointed products) that overlap in time and in the sky, as seen
    from Solar Orbiter: columns block_a, block_b, start, end (of the common time), frac
    (fraction of block_a's FOV inside block_b's)."""
    from matplotlib.path import Path as MplPath
    out = []
    for _, a in blocks_a.iterrows():
        common = blocks_b[(blocks_b.start <= a.end) & (blocks_b.end >= a.start)]
        for _, b in common.iterrows():
            frac = np.nan
            if "c0x" in a and "c0x" in b and np.isfinite(a.c0x) and np.isfinite(b.c0x):
                ca = a[CORNER_COLS].to_numpy(float).reshape(4, 2)
                cb = b[CORNER_COLS].to_numpy(float).reshape(4, 2)
                frac = MplPath(cb).contains_points(fov_grid(ca, 30)).mean()
            if np.isnan(frac) or frac > min_frac:
                out.append(dict(block_a=a.block, block_b=b.block, start=max(a.start, b.start),
                                end=min(a.end, b.end), frac=frac))
    return pd.DataFrame(out)


def download(filename, folder, overwrite=False):
    """Download one Solar Orbiter file into folder (kept: not downloaded again); returns the
    path. EUI files come from the EUI team's server at ROB (data.observatory.be: the same
    files, ~6x faster than SOAR on 2025-11-04), anything else or anything missing there
    from SOAR (latest version of the data item)."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / filename
    if path.is_file() and not overwrite:
        return path
    r = None
    m = re.match(r"solo_(L\d)_eui-.*?_(\d{4})(\d{2})(\d{2})T", filename)
    if m:
        url = f"{EUI_URL}/{m.group(1)}/{m.group(2)}/{m.group(3)}/{m.group(4)}/{filename}"
        try:
            r = requests.get(url, timeout=600, stream=True)
            if r.status_code != 200:
                r = None
        except requests.RequestException:
            r = None
    if r is None:
        did = re.sub(r"_V\d+.*$", "", filename)
        r = requests.get(DATA_URL, params=dict(retrieval_type="LAST_PRODUCT", data_item_id=did,
                                               product_type="SCIENCE"), timeout=600, stream=True)
    r.raise_for_status()
    tmp = path.with_suffix(".part")
    with open(tmp, "wb") as f:
        for chunk in r.iter_content(1 << 20):
            f.write(chunk)
    tmp.replace(path)
    return path
