#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
Where Solar Orbiter's high-resolution telescopes (HRIEUV, HRI-Lya, PHI-HRT) and DKIST point
at a given time, where PSP (and other spacecraft) are, and which part of the Sun their
solar wind can come from, drawn on the full disk from both vantage points:

  left    AIA (Earth view; 1024x1024 synoptic image, fetched over HTTP in a second)
  right   EUI/FSI (Solar Orbiter's view)
  (third  with --positions: the spacecraft seen from above the ecliptic, as in
          plot_positions.py, with the Parker spirals of the --connect spacecraft)

Solar Orbiter FOVs
  on FSI  as observed: solid on the disk, dashed off the limb
  on AIA  the part on Solar Orbiter's disk as its footprint on the solar surface (solid where
          AIA sees that surface, dotted where it is behind the disk, Earth's far side). The
          part off Solar Orbiter's limb sees the corona: for HRIEUV, HRI-Lya and SPICE it is
          drawn dashed where that corona appears from Earth, on a sphere through the Sun's
          centre seen from Solar Orbiter (a guide: the depth of the emission along the line of
          sight is unknown); PHI-HRT sees no corona, so its outline follows Solar Orbiter's limb.
DKIST FOVs (DKIST Data Center, dkist_meta.py: bounding box of each dataset overlapping
  [time - --dkist-window, time + --dkist-window], Earth view): as observed on AIA; on FSI
  the same way as the coronal Solar Orbiter FOVs on AIA (footprint, and the off-limb part
  on the screen, dashed).
Spacecraft (--spacecraft, default PSP; positions from JPL Horizons, shared cache with
  plot_positions.py): sub-spacecraft point (x; o at the limb when on the far side) and the
  spacecraft's limb, on both images. Solar Orbiter is always shown on AIA, Earth on FSI.
Connectivity region (--connect, default PSP), as in plot_ref_file.py: the spacecraft's
  latitude, and the longitudes a Parker spiral for --speeds (default 300-800 km/s) reaches
  at 2.5 Rsun west of the sub-spacecraft point, widened by --margin (default 10 deg); the
  line along its middle is labelled with the speeds at its ends. Ballistic, no PFSS: in a
  coronal hole the photospheric footpoint can lie elsewhere.

Pointing from metadata only (soar_meta.py, dkist_meta.py): only the background images are
downloaded (AIA ~1 MB, FSI ~4 MB, cached in cache/). For each Solar Orbiter product the
file nearest in time (within --window) is used.

Run:  ./plot_fov.py 2025-10-20T12:30
      ./plot_fov.py 2025-03-16T12:00 --aia 171 --fsi 304 --zoom
      ./plot_fov.py 2025-10-30T19:30 --window 12h --positions     # + DKIST Cryo-NIRSP
      ./plot_fov.py 2025-10-20T12:30 --connect PSP "Solar Orbiter" "Earth / L1" --speeds 350-750
      ./plot_fov.py 2025-10-20T12:30 --no-show            # only save fov/fov_20251020_1230.png
"""
import argparse
import textwrap
import warnings
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval
from matplotlib import patheffects
from matplotlib.lines import Line2D
from sunpy.coordinates import HeliographicStonyhurst, get_earth
from sunpy.coordinates.sun import L0 as carrington_L0
from sunpy.visualization.drawing import limb

import dkist_meta as dm
import soar_meta as sm
from plot_positions import R_SS, SPACECRAFT, V_RANGE, Trajectories, parse_speeds, spiral_shift

OUT_DIR = Path(__file__).resolve().parent
FIG_DIR = OUT_DIR / "fov"
AIA_CACHE = OUT_DIR / "cache" / "aia_synoptic"
FSI_CACHE = OUT_DIR / "cache" / "fsi"
AIA_URL = "http://jsoc1.stanford.edu/data/aia/synoptic/{t:%Y/%m/%d}/H{t:%H}00/AIA{t:%Y%m%d_%H%M}_{wl:04d}.fits"
# one colour per instrument (both HRIEUV products, 174 and filter wheel "open", share it)
COLORS = {"hri174": "lime", "hriopn": "lime", "hrilya": "hotpink", "phi-hrt": "cyan", "spice": "gold"}
LABELS = {"hri174": "HRIEUV", "hriopn": "HRIEUV", "hrilya": "HRI-Lya", "phi-hrt": "PHI-HRT", "spice": "SPICE"}
# spacecraft as in plot_ref_file.py; Earth / L1 (Wind, ACE): connectivity region only
SC_COLORS = {"PSP": "magenta", "Solar Orbiter": "orange", "STEREO-A": "yellow", "BepiColombo": "dodgerblue",
             "SOHO": "white", "SDO": "white", "Earth / L1": "red"}
L1 = "Earth / L1"
EARTH_COLOR = "deepskyblue"
DKIST_COLORS = {"VBI": "white", "VISP": "red", "DL-NIRSP": "cornflowerblue", "CRYO-NIRSP": "mediumorchid", "VTF": "chocolate"}
DKIST_OTHER = "silver"
IRIS_COLORS = {"raster": "coral", "sji": "peachpuff"}    # spectrograph raster area, slit-jaw area (browse_web)
assert len(set(COLORS.values()) | set(DKIST_COLORS.values()) | {DKIST_OTHER} | set(IRIS_COLORS.values())) == \
    len(set(LABELS.values())) + len(DKIST_COLORS) + 1 + len(IRIS_COLORS), "two instruments share a FOV colour"
OUTLINE = [patheffects.withStroke(linewidth=2.5, foreground="black")]
N_EDGE = 80     # points per FOV edge, so the reprojected edges can bend
MARGIN = 10     # deg, default margin around a connectivity region


# ---------------------------------------------------------------- background images

def aia_synoptic(t, wavelength):
    """The AIA synoptic file nearest to t (2-min cadence), downloaded to cache/aia_synoptic/."""
    AIA_CACHE.mkdir(parents=True, exist_ok=True)
    t0 = pd.Timestamp(t).round("2min")
    for dt in (0, 2, -2, 4, -4, 6, -6, 10, -10, 20, -20):
        tt = t0 + pd.Timedelta(minutes=dt)
        url = AIA_URL.format(t=tt, wl=wavelength)
        path = AIA_CACHE / Path(url).name
        if path.is_file():
            return path
        r = requests.get(url, timeout=120)
        if r.status_code == 200 and r.content[:6] == b"SIMPLE":
            path.write_bytes(r.content)
            return path
    raise FileNotFoundError(f"no AIA {wavelength} synoptic file within 20 min of {t0}")


def fsi_nearest(t, band, windows=("1h", "6h", "1D")):
    """The FSI file (band 174 / 304) nearest to t, looking within each window in turn,
    downloaded to cache/fsi/."""
    t = pd.Timestamp(t)
    for window in windows:
        w = pd.Timedelta(window)
        f = sm.files(t - w, t + w, [f"fsi{band}"], headers=False)
        if len(f):
            row = f.iloc[np.argmin(np.abs((f.start - t).dt.total_seconds()))]
            if abs(row.start - t) > pd.Timedelta("10min"):
                print(f"FSI {band}: nearest file {row.start:%Y-%m-%d %H:%M} ({(row.start - t).total_seconds() / 3600:+.1f} h)")
            return sm.download(row.filename, FSI_CACHE)
    return None


def show_map(ax, m, extent=1.25):
    """Map with its own colour map, asinh stretch; limits +-extent solar radii."""
    data = np.where(np.isfinite(m.data), m.data, 0)
    norm = ImageNormalize(data, interval=PercentileInterval(99.7), stretch=AsinhStretch(0.05))
    m.plot(axes=ax, norm=norm, title=False)
    m.draw_grid(axes=ax, grid_spacing=15 * u.deg, color="white", alpha=0.25, lw=0.5)
    ax.set_facecolor("black")
    r = extent * m.rsun_obs.to_value(u.arcsec)
    corners = SkyCoord([-r, r] * u.arcsec, [-r, r] * u.arcsec, frame=m.coordinate_frame)
    x, y = m.wcs.world_to_pixel(corners)
    ax.set_xlim(*x)
    ax.set_ylim(*y)


# ---------------------------------------------------------------- geometry and drawing

def box_outline(corners):
    """Densely sampled closed outline of a quadrilateral (4, 2) -> tx, ty."""
    c = np.vstack([corners, corners[:1]])
    s = np.linspace(0, 1, N_EDGE, endpoint=False)
    pts = np.vstack([c[k] + s[:, None] * (c[k + 1] - c[k]) for k in range(4)] + [c[:1]])
    return pts[:, 0], pts[:, 1]


def to_surface(tx, ty, lat_obs, lon_obs, r_au):
    """Outline (helioprojective arcsec seen from an observer at Stonyhurst lat_obs, lon_obs, r_au)
    -> Stonyhurst lon, lat of its footprint on the surface: points off the observer's limb
    are moved onto the limb, in their direction."""
    tx, ty = np.asarray(tx, float), np.asarray(ty, float)
    alpha = np.degrees(np.arcsin(sm.RSUN_KM / (r_au * sm.AU_KM))) * 3600   # limb radius, arcsec
    rho = np.hypot(tx, ty)
    f = np.where(rho > alpha, alpha / np.where(rho > 0, rho, 1), 1.0)
    lon, lat = sm.hpc_to_hgs(tx * f, ty * f, lat_obs, lon_obs, r_au)
    shrink = 1.0
    while np.isnan(lon).any() and shrink > 0.99:   # just past the limb by rounding: step inside
        shrink -= 1e-4
        bad = np.isnan(lon)
        lon[bad], lat[bad] = sm.hpc_to_hgs(tx[bad] * f[bad] * shrink, ty[bad] * f[bad] * shrink,
                                           lat_obs, lon_obs, r_au)
    return lon, lat


def hgs_to_hpc(lon, lat, lat_obs, lon_obs, r_au, r_km=None):
    """Stonyhurst lon, lat (deg) at distance r_km from the Sun's centre (default: on the
    surface) -> helioprojective tx, ty (arcsec) seen from an observer at Stonyhurst lat_obs,
    lon_obs, r_au, and z (km, the heliocentric coordinate towards the observer; < 0 behind the
    plane of the sky). Closed form, the inverse of soar_meta.hpc_to_hgs."""
    D, R = r_au * sm.AU_KM, sm.RSUN_KM if r_km is None else np.asarray(r_km, float)
    B0, L0 = np.radians(lat_obs), np.radians(lon_obs)
    lon, lat = np.radians(np.asarray(lon, float)), np.radians(np.asarray(lat, float))
    x = R * np.cos(lat) * np.sin(lon - L0)
    y = R * (np.sin(lat) * np.cos(B0) - np.cos(lat) * np.cos(lon - L0) * np.sin(B0))
    z = R * (np.sin(lat) * np.sin(B0) + np.cos(lat) * np.cos(lon - L0) * np.cos(B0))
    d = np.sqrt(x ** 2 + y ** 2 + (D - z) ** 2)
    return np.degrees(np.arctan2(x, D - z)) * 3600, np.degrees(np.arcsin(y / d)) * 3600, z


def hpc_to_screen(tx, ty, lat_obs, lon_obs, r_au):
    """Lines of sight (helioprojective arcsec) of an observer -> Stonyhurst lon, lat (deg)
    and distance from the Sun's centre (km) of the point where each meets a sphere centred
    on the observer through the Sun's centre (sunpy's SphericalScreen): where off-limb
    emission is put when its depth along the line of sight is unknown."""
    D = r_au * sm.AU_KM
    B0, L0 = np.radians(lat_obs), np.radians(lon_obs)
    tx, ty = np.radians(np.asarray(tx, float) / 3600), np.radians(np.asarray(ty, float) / 3600)
    x, y, z = D * np.cos(ty) * np.sin(tx), D * np.sin(ty), D - D * np.cos(ty) * np.cos(tx)
    xp, zz = z * np.cos(B0) - y * np.sin(B0), z * np.sin(B0) + y * np.cos(B0)
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    return np.degrees(L0 + np.arctan2(x, xp)), np.degrees(np.arcsin(zz / r)), r


# Solar Orbiter products that see the corona above the limb (their off-limb part is drawn going up off
# the limb on the other view); PHI-HRT sees the photosphere only (off-limb pixels hold no data: footprint only)
CORONAL = {"hri174", "hriopn", "hrilya", "spice"}


def project(tx, ty, src, dst, coronal=True):
    """An outline observed from src = (lat, lon, r_au) (Stonyhurst deg, au) as seen from dst:
    x, y (arcsec in dst's view) and kind per point:
      on      on src's disk: its footprint on the surface, which dst sees
      hidden  on src's disk: its footprint on the surface, on dst's far side
      off     off src's limb (corona): on the screen through the Sun's centre (hpc_to_screen);
              'hidden' instead where the Sun hides it from dst
    coronal=False (photosphere only): the off-limb part is moved onto src's limb instead,
    so the footprint follows src's limb (kinds on / hidden only)."""
    tx, ty = np.asarray(tx, float), np.asarray(ty, float)
    if not coronal:
        lon, lat = to_surface(tx, ty, *src)
        ok = np.isfinite(lon) & np.isfinite(lat)
        x, y, z = hgs_to_hpc(lon[ok], lat[ok], *dst)
        vis = cos_seen(lon[ok], lat[ok], dst)
        return x, y, np.where(vis, "on", "hidden")
    lon, lat = sm.hpc_to_hgs(tx, ty, *src)
    on = np.isfinite(lon)
    x, y = np.full(len(tx), np.nan), np.full(len(tx), np.nan)
    kind = np.full(len(tx), "off", dtype=object)
    if on.any():
        x[on], y[on], _ = hgs_to_hpc(lon[on], lat[on], *dst)
        kind[on] = np.where(cos_seen(lon[on], lat[on], dst), "on", "hidden")
    if (~on).any():
        slon, slat, sr = hpc_to_screen(tx[~on], ty[~on], *src)
        xo, yo, zo = hgs_to_hpc(slon, slat, *dst, r_km=sr)
        x[~on], y[~on] = xo, yo
        alpha = np.degrees(np.arcsin(sm.RSUN_KM / (dst[2] * sm.AU_KM))) * 3600   # dst's limb radius
        kind[~on] = np.where((zo < 0) & (np.hypot(xo, yo) < alpha), "hidden", "off")
    ok = np.isfinite(x) & np.isfinite(y)
    return x[ok], y[ok], kind[ok].astype(str)


def cos_seen(lon, lat, obs):
    """True where surface points face an observer obs = (lat, lon, r_au)."""
    return sm.cos_angle(lon, lat, obs[1], obs[0]) > sm.RSUN_KM / (obs[2] * sm.AU_KM)


def visible_from(m, lon, lat):
    """True where Stonyhurst lon, lat on the surface face the observer of map m."""
    o = m.observer_coordinate
    return sm.cos_angle(lon, lat, o.lon.to_value(u.deg), o.lat.to_value(u.deg)) > \
        (m.rsun_meters / m.dsun).decompose().value


def draw_runs(ax, x, y, kind, color, styles, lw=1.8, **kw):
    """Plot a polyline in runs of equal kind, each with styles[kind] (None = skip)."""
    start = 0
    for k in range(1, len(kind) + 1):
        if k == len(kind) or kind[k] != kind[start]:
            ls = styles.get(kind[start])
            if ls is not None:
                sl = slice(max(start - 1, 0), k)   # overlap by one point: no gaps
                ax.plot(x[sl], y[sl], ls, color=color, lw=lw, path_effects=OUTLINE, **kw)
            start = k


def draw_surface(ax, m, lon, lat, color, lw=1.8):
    """A line on the solar surface (Stonyhurst lon, lat) on map m: solid where m's observer
    sees it, dotted where it is behind the disk. Returns the pixel coordinates."""
    ok = np.isfinite(lon) & np.isfinite(lat)
    lon, lat = lon[ok], lat[ok]
    c = SkyCoord(lon * u.deg, lat * u.deg, m.rsun_meters, frame=HeliographicStonyhurst(obstime=m.date))
    px = m.world_to_pixel(c)
    x, y = px.x.value, px.y.value
    kind = np.where(visible_from(m, lon, lat), "on", "hidden")
    draw_runs(ax, x, y, kind, color, {"on": "-", "hidden": ":"}, lw=lw)
    return x, y


def draw_hpc_box(ax, m, tx, ty, lat_obs, lon_obs, r_au, color, lw=1.8):
    """An outline observed from the map's own vantage point: solid on the disk, dashed off it."""
    px = m.wcs.world_to_pixel(SkyCoord(tx * u.arcsec, ty * u.arcsec, frame=m.coordinate_frame))
    lon, _ = sm.hpc_to_hgs(tx, ty, lat_obs, lon_obs, r_au)
    draw_runs(ax, px[0], px[1], np.where(np.isfinite(lon), "on", "off"), color, {"on": "-", "off": "--"}, lw=lw)
    return px


def map_observer(m):
    """Stonyhurst lat, lon (deg) and distance (au) of the observer of map m."""
    o = m.observer_coordinate
    return o.lat.to_value(u.deg), o.lon.to_value(u.deg), o.radius.to_value(u.au)


# ---------------------------------------------------------------- Solar Orbiter FOVs

def nearest_files(t, products, window):
    """For each product, the file nearest to t within window (a row of soar_meta.files())."""
    t, w = pd.Timestamp(t), pd.Timedelta(window)
    f = sm.files(t - w, t + w, products)
    out = {}
    for product in products:
        g = f[(f["product"] == product)].dropna(subset=["c0x"])
        if len(g):
            out[product] = g.iloc[np.argmin(np.abs((g.start - t).dt.total_seconds()))]
    return out


def solo_observer(row):
    return SkyCoord(row.solo_lon * u.deg, row.solo_lat * u.deg, row.solo_r * u.au,
                    frame=HeliographicStonyhurst(obstime=row.start))


def fov_outline(row):
    return box_outline(row[sm.CORNER_COLS].to_numpy(float).reshape(4, 2))


def draw_fov_solo(ax, m, row, color):
    """FOV as observed, on a Solar Orbiter (FSI) map."""
    tx, ty = fov_outline(row)
    return draw_hpc_box(ax, m, tx, ty, row.solo_lat, row.solo_lon, row.solo_r, color)


PROJECTED = {"on": "-", "hidden": ":", "off": "--"}   # line styles of project()'s kinds


def draw_projected(ax, m, tx, ty, src, color, coronal, lw=1.8):
    """An outline observed from src = (lat, lon, r_au) drawn on map m (another vantage point),
    see project(): solid on the surface m sees, dotted hidden, dashed off src's limb.
    Returns the pixel coordinates."""
    x, y, kind = project(tx, ty, src, map_observer(m), coronal)
    if not len(x):
        return np.array([]), np.array([])
    px = m.wcs.world_to_pixel(SkyCoord(x * u.arcsec, y * u.arcsec, frame=m.coordinate_frame))
    draw_runs(ax, px[0], px[1], kind, color, PROJECTED, lw=lw)
    return px


def draw_fov_earth(ax, m, row, color):
    """A Solar Orbiter FOV on the AIA map (see the module docstring)."""
    tx, ty = fov_outline(row)
    return draw_projected(ax, m, tx, ty, (row.solo_lat, row.solo_lon, row.solo_r), color,
                          row["product"] in CORONAL)


# ---------------------------------------------------------------- DKIST

def dkist_fovs(t, window):
    """DKIST datasets overlapping [t - window, t + window], one per instrument and box."""
    t, w = pd.Timestamp(t), pd.Timedelta(window)
    try:
        d = dm.datasets(t - w, t + w)
    except requests.RequestException as e:
        print(f"DKIST: {e}")
        return pd.DataFrame()
    if d.empty:
        return d
    d = d.dropna(subset=["x0"])
    key = d[["x0", "y0", "x1", "y1"]].round(-1).astype(str).agg(",".join, axis=1)
    return d.assign(key=d.instrument + key).drop_duplicates("key").reset_index(drop=True)


def draw_dkist(ax_aia, aia, ax_fsi, fsi, d):
    earth = map_observer(aia)
    for _, r in d.iterrows():
        color = DKIST_COLORS.get(str(r.instrument).upper(), DKIST_OTHER)
        corners = np.array([[r.x0, r.y0], [r.x1, r.y0], [r.x1, r.y1], [r.x0, r.y1]])
        tx, ty = box_outline(corners)
        draw_hpc_box(ax_aia, aia, tx, ty, *earth, color, lw=1.4)
        if ax_fsi is not None:   # DKIST sees the corona off the limb too (e.g. Cryo-NIRSP)
            draw_projected(ax_fsi, fsi, tx, ty, earth, color, coronal=True, lw=1.4)


# ---------------------------------------------------------------- spacecraft

def positions(t, names, rows):
    """Stonyhurst coordinates of the spacecraft names (and Earth / L1) at t."""
    t = Time(t)
    traj = Trajectories()
    traj.fetch(t, [n for n in names if n in SPACECRAFT], status=print)
    out = {}
    for n in names:
        if n == L1:
            e = get_earth(t)
            out[n] = SkyCoord(e.lon, e.lat, 0.99 * e.radius, frame=HeliographicStonyhurst(obstime=t))
        else:
            c = traj.position(n, t)
            if c is None and n == "Solar Orbiter" and rows:
                c = solo_observer(next(iter(rows.values())))
            out[n] = None if c is None else SkyCoord(c)
    return out


def mark_observer(ax, m, obs, color, name):
    """Sub-observer point of obs on map m (open circle at the limb if on the far side) and
    obs's limb (solid where m's observer sees it, dotted where hidden)."""
    sub = SkyCoord(obs.lon, obs.lat, m.rsun_meters, frame=HeliographicStonyhurst(obstime=m.date))
    front = visible_from(m, obs.lon.to_value(u.deg), obs.lat.to_value(u.deg))
    hpc = sub.transform_to(m.coordinate_frame)
    if not front:   # at the limb, in the direction of the hidden sub-point
        a = np.arctan2(hpc.Ty.to_value(u.arcsec), hpc.Tx.to_value(u.arcsec))
        r = m.rsun_obs.to_value(u.arcsec)
        hpc = SkyCoord(r * np.cos(a) * u.arcsec, r * np.sin(a) * u.arcsec, frame=m.coordinate_frame)
    x, y = m.wcs.world_to_pixel(hpc)
    ax.plot(x, y, "x" if front else "o", color=color, ms=9, mew=2, mfc="none", path_effects=OUTLINE)
    ax.annotate(name if front else f"{name} (far side)", (x, y), xytext=(6, 5), textcoords="offset points",
                color=color, fontsize=9, path_effects=OUTLINE)
    limb(ax, SkyCoord(obs), rsun=m.rsun_meters, color=color, lw=1.2)


def connectivity(c, v_range, margin):
    """Connectivity region of a spacecraft at c (see the module docstring): Stonyhurst
    longitudes reached at R_SS for the fast and slow speed, and the region's outline lon, lat."""
    lon0, lat0 = c.lon.to_value(u.deg), c.lat.to_value(u.deg)
    fast, slow = (lon0 + spiral_shift(c.radius, R_SS, v) for v in v_range[::-1])
    lon_a, lon_b = fast - margin, slow + margin
    lat_a, lat_b = max(lat0 - margin, -89.9), min(lat0 + margin, 89.9)
    n = 60   # points per edge, so the edges follow the sphere
    lon = np.r_[np.linspace(lon_a, lon_b, n), np.full(n, lon_b), np.linspace(lon_b, lon_a, n), np.full(n, lon_a)]
    lat = np.r_[np.full(n, lat_a), np.linspace(lat_a, lat_b, n), np.full(n, lat_b), np.linspace(lat_b, lat_a, n)]
    return fast, slow, lon, lat


def draw_connectivity(ax, m, c, name, color, v_range, margin):
    fast, slow, lon, lat = connectivity(c, v_range, margin)
    x, y = draw_surface(ax, m, lon, lat, color, lw=2)
    seen = visible_from(m, lon, lat)
    top = 2 * len(lon) // 4 + len(lon) // 8   # middle of the north edge
    if seen[top]:
        ax.annotate(name, (x[top], y[top]), xytext=(0, 4), textcoords="offset points", ha="center", va="bottom",
                    color=color, fontsize=9, path_effects=OUTLINE)
    # the speed range itself: a line along the spacecraft's latitude, ends labelled in km/s
    line_lon = np.linspace(fast, slow, 40)
    line_lat = np.full_like(line_lon, c.lat.to_value(u.deg))
    x, y = draw_surface(ax, m, line_lon, line_lat, color, lw=2.5)
    ends = visible_from(m, line_lon[[0, -1]], line_lat[[0, -1]])
    for k, v, ha in ((0, v_range[1], "right"), (-1, v_range[0], "left")):
        if ends[k]:
            ax.annotate(f"{v:g}", (x[k], y[k]), xytext=(0, -9), textcoords="offset points", ha=ha,
                        color=color, fontsize=8, path_effects=OUTLINE)
    return fast, slow


def draw_positions(ax, t, coords, connect, v_range):
    """Polar plot (Stonyhurst, seen from above the ecliptic; Earth at lon 0) of the
    spacecraft, with the Parker spirals of the connect ones (as plot_positions.py)."""
    earth = get_earth(Time(t))
    objs = {"Earth": SkyCoord(earth)} | {n: c for n, c in coords.items() if c is not None and n != L1}
    ax.plot(0, 0, "o", color="orange", ms=10)
    rmax = 1.1 * max(c.radius.to_value(u.au) for c in objs.values())
    for n, c in objs.items():
        color = EARTH_COLOR if n == "Earth" else SC_COLORS.get(n, "k")
        lon, r = c.lon.to_value(u.rad), c.radius.to_value(u.au)
        if n in connect or (n == "Earth" and L1 in connect):
            rr = np.linspace(R_SS.to_value(u.au), r, 200) * u.au
            th = [np.radians(c.lon.to_value(u.deg) + spiral_shift(c.radius, rr, v)) for v in v_range]
            ax.fill(np.r_[th[0], th[1][::-1]], np.r_[rr.value, rr.value[::-1]], color=color, alpha=0.15, lw=0)
            for t_, ls in zip(th, ("--", ":")):
                ax.plot(t_, rr.value, ls, color=color, lw=0.9)
        ax.plot([0, lon], [0, r], color=color, lw=0.8, alpha=0.6)
        ax.plot(lon, r, "o" if n == "Earth" else "D", color=color, ms=8 if n == "Earth" else 6, mec="k", mew=0.5)
        ax.annotate(n, (lon, r), xytext=(5, 4), textcoords="offset points", fontsize=8, color="k")
    ax.set_rmax(rmax)
    ax.set_rlabel_position(-22.5)
    ax.set_thetagrids(range(0, 360, 45), [f"{a if a <= 180 else a - 360}°" for a in range(0, 360, 45)])
    ax.tick_params(labelsize=8)
    ax.set_title("positions (Stonyhurst, au; Earth at 0°)", fontsize=10, pad=14)


# ---------------------------------------------------------------- figure

def legend(ax, items):
    """items: [(label, color, linestyle)]"""
    handles = [Line2D([], [], color=c, ls=ls, lw=2, path_effects=OUTLINE, label=lab) for lab, c, ls in items]
    if handles:
        leg = ax.legend(handles=handles, loc="lower left", fontsize=8, framealpha=0.6, facecolor="black")
        for text in leg.get_texts():
            text.set_color("white")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("time", help="UTC, e.g. 2025-10-20T12:30")
    ap.add_argument("--products", nargs="*", default=["hri174", "phi-hrt"], choices=sm.POINTED)
    ap.add_argument("--window", default="2h", help="Solar Orbiter files within this of the time (default 2h)")
    ap.add_argument("--aia", type=int, default=193, help="AIA channel (default 193)")
    ap.add_argument("--fsi", type=int, default=174, choices=[174, 304], help="FSI channel (default 174)")
    ap.add_argument("--no-fsi", action="store_true", help="AIA panel only")
    ap.add_argument("--spacecraft", nargs="*", default=["PSP"], choices=[n for n in SPACECRAFT if n != "Solar Orbiter"],
                    help="mark these (sub-spacecraft point + limb; default PSP)")
    ap.add_argument("--connect", nargs="*", default=["PSP"], choices=list(SC_COLORS),
                    help='connectivity regions of these (default PSP; e.g. PSP "Solar Orbiter" "Earth / L1")')
    ap.add_argument("--speeds", default=f"{V_RANGE[0]}-{V_RANGE[1]}", help="wind speed range, km/s (default 300-800)")
    ap.add_argument("--margin", type=float, default=MARGIN, help="deg around the connectivity region (default 10)")
    ap.add_argument("--no-dkist", action="store_true", help="do not look for DKIST datasets")
    ap.add_argument("--dkist-window", default="3h", help="DKIST datasets within this of the time (default 3h)")
    ap.add_argument("--positions", action="store_true", help="add a panel with the spacecraft positions")
    ap.add_argument("--zoom", action="store_true", help="zoom both images on the Solar Orbiter FOVs")
    ap.add_argument("--no-show", action="store_true", help="only save the figure")
    args = ap.parse_args()
    t = pd.Timestamp(args.time)
    v_range = parse_speeds(args.speeds)
    warnings.filterwarnings("ignore", module="sunpy")
    warnings.filterwarnings("ignore", module="astropy")

    rows = nearest_files(t, args.products, args.window) if args.products else {}
    for p in args.products:
        if p in rows:
            r = rows[p]
            print(f"{LABELS[p]:8} {r.filename}  ({r.soop}, pointing from {r.pointing_src})")
        else:
            print(f"{LABELS[p]:8} no file within {args.window}")
    dk = pd.DataFrame() if args.no_dkist else dkist_fovs(t, args.dkist_window)
    for _, r in dk.iterrows():
        print(f"DKIST    {r.instrument:10} {r.start:%m-%d %H:%M}-{r.end:%H:%M}  {r.proposal}  {r.lines}")
    names = list(dict.fromkeys(["Solar Orbiter"] + args.spacecraft + args.connect))
    coords = positions(t, names, rows)

    aia = sunpy.map.Map(aia_synoptic(t, args.aia))
    fsi_path = None if args.no_fsi else fsi_nearest(t, args.fsi)
    fsi = sunpy.map.Map(fsi_path) if fsi_path else None
    if fsi_path is None and not args.no_fsi:
        print(f"FSI {args.fsi}: no file within 1 day")

    widths = [1] + ([1] if fsi is not None else []) + ([0.75] if args.positions else [])
    fig = plt.figure(figsize=(7.2 * sum(widths), 7.8))
    gs = fig.add_gridspec(1, len(widths), width_ratios=widths)
    ax1 = fig.add_subplot(gs[0], projection=aia)
    show_map(ax1, aia)
    ax1.set_title(f"SDO/AIA {args.aia} Å  {aia.date.strftime('%Y-%m-%d %H:%M:%S')}  (Earth view)", fontsize=10)
    ax2 = None
    if fsi is not None:
        ax2 = fig.add_subplot(gs[1], projection=fsi)
        show_map(ax2, fsi)
        ax2.set_title(f"SolO/EUI FSI {args.fsi} Å  {fsi.date.strftime('%Y-%m-%d %H:%M:%S')}  "
                      f"(Solar Orbiter view)", fontsize=10)

    items, boxes1, boxes2 = [], [], []
    for p, r in rows.items():
        boxes1.append(draw_fov_earth(ax1, aia, r, COLORS[p]))
        if ax2 is not None:
            boxes2.append(draw_fov_solo(ax2, fsi, r, COLORS[p]))
        items.append((LABELS[p], COLORS[p], "-"))
    if len(dk):
        draw_dkist(ax1, aia, ax2, fsi, dk)
        items += [(f"DKIST {i}", DKIST_COLORS.get(str(i).upper(), DKIST_OTHER), "-") for i in dict.fromkeys(dk.instrument)]

    info = []
    carrington = carrington_L0(Time(t)).to_value(u.deg)
    earth = get_earth(aia.date)
    for n in names:
        c = coords.get(n)
        if c is None:
            info.append(f"{n}: no position")
            continue
        color = SC_COLORS.get(n, "white")
        line = f"{n}:"
        if n != L1:
            sep = np.degrees(np.arccos(np.clip(sm.cos_angle(c.lon.to_value(u.deg), c.lat.to_value(u.deg),
                                                            0.0, earth.lat.to_value(u.deg)), -1, 1)))
            line += (f" r {c.radius.to_value(u.au):.3f} au ({c.radius.to_value(u.R_sun):.0f} Rsun), "
                     f"lon {c.lon.wrap_at(180 * u.deg).to_value(u.deg):+.1f}, lat {c.lat.to_value(u.deg):+.1f}, "
                     f"{sep:.0f}° from Earth")
            if n == "Solar Orbiter" or n in args.spacecraft:
                mark_observer(ax1, aia, c, color, n)
                if ax2 is not None and n != "Solar Orbiter":
                    mark_observer(ax2, fsi, c, color, n)
        if n in args.connect:
            for ax, m in ((ax1, aia), (ax2, fsi)):
                if ax is not None:
                    fast, slow = draw_connectivity(ax, m, c, n, color, v_range, args.margin)
            line += (f"; connects lon {fast:+.0f}..{slow:+.0f} (Carr. {(fast + carrington) % 360:.0f}.."
                     f"{(slow + carrington) % 360:.0f}) ±{args.margin:g}° for {v_range[0]:g}-{v_range[1]:g} km/s")
            items.append((f"{n} source region", color, "-"))
        info.append(line)
    if ax2 is not None:
        mark_observer(ax2, fsi, earth, EARTH_COLOR, "Earth")
    legend(ax1, items)
    if ax2 is not None:
        legend(ax2, items)
    if args.positions:
        draw_positions(fig.add_subplot(gs[-1], projection="polar"), t, coords, args.connect, v_range)

    if args.zoom and boxes1:
        for ax, boxes in ((ax1, boxes1), (ax2, boxes2)):
            if ax is None or not boxes:
                continue
            x = np.concatenate([b[0] for b in boxes]); y = np.concatenate([b[1] for b in boxes])
            ok = np.isfinite(x) & np.isfinite(y)
            pad = 0.15 * max(np.ptp(x[ok]), np.ptp(y[ok]))
            ax.set_xlim(x[ok].min() - pad, x[ok].max() + pad)
            ax.set_ylim(y[ok].min() - pad, y[ok].max() + pad)

    fig.suptitle(f"{t:%Y-%m-%d %H:%M} UTC", fontsize=11)
    lines = [f"{LABELS[p]}: {r.start:%H:%M:%S}  {sm.PRODUCTS[p]['desc'][0]}  {str(r.soop).split('_')[-1]}"
             for p, r in rows.items()] + info
    lines.append("Solar Orbiter FOVs on AIA: footprint on the surface (solid seen, dotted behind the disk), "
                 "off-limb part (corona) dashed on a sphere through the Sun's centre (PHI-HRT: footprint only); "
                 "on FSI: as observed (dashed off limb). Source regions: ballistic Parker spiral to 2.5 Rsun, no PFSS")
    width = int(fig.get_figwidth() * 15.5)   # characters per line at this font size
    lines = [w for line in lines for w in textwrap.wrap(line, width, subsequent_indent="    ")]
    fig.text(0.01, 0.01, "\n".join(lines), fontsize=7.5, family="monospace", va="bottom")
    fig.tight_layout(rect=(0, 0.03 + 0.016 * len(lines), 1, 0.97))

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / f"fov_{t:%Y%m%d_%H%M}.png"
    fig.savefig(out, dpi=130)
    print(f"Saved {out}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
