#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
The time browser in a web page (Panel + Bokeh on localhost): what was observed when, and
where. It replaces browse.py (same data, same cache), faster and with zoom / pan / hover.

  timeline   one row per instrument (HRIEUV, HRI-Lya, PHI-HRT, SPICE, FSI, DKIST): a bar per
             observing block (Solar Orbiter), FSI run or DKIST dataset of the months loaded so
             far; hover for details, click to go to that time; the shaded band is the step.
             Wheel zooms the time axis, drag pans it.
  AIA        (Earth view; synoptic 1024x1024, 2-min cadence) with the Solar Orbiter FOVs:
             the part on Solar Orbiter's disk as its footprint on the surface (solid where AIA
             sees it, dotted behind the disk); the part off Solar Orbiter's limb (corona;
             HRIEUV, HRI-Lya, SPICE) dashed, going up off the limb (plot_fov.project; PHI-HRT
             sees no corona: its outline follows Solar Orbiter's limb); DKIST FOVs, Solar Orbiter's sub-point and limb, and a 15° grid
  FSI 174    (Solar Orbiter view, nearest image within a day; shown when HRI / PHI-HRT /
             SPICE observed in the step) with the FOVs as observed (dashed off the limb), DKIST
             FOVs the same way as the coronal ones on AIA, Earth's sub-point and limb
             Both images: wheel zooms, drag pans, box zoom, reset; the zoom is kept when
             stepping. Below each image: the cursor position, helioprojective x, y (arcsec) and
             Stonyhurst lon, lat and mu (off the disk: the distance from the centre in Rsun).
             Hover inside a FOV for its block: time range, cadence, files, size, where.
  tabs       summary (what each instrument measures, cadence, time range, SOOP and its aim,
             SPICE study and lines, FSI, DKIST; as in browse.py), Solar Orbiter files of the
             step, DKIST datasets of the step
  controls   date (UTC, Enter), << / >> by the cadence, cadence, AIA channel, which
             observations to draw, and
    positions                    (on / off) PSP and Solar Orbiter on both images: sub-point (x;
                                 o at the limb when on the far side), limb, and solar source
                                 region (ballistic Parker spiral to 2.5 Rsun, 300-800 km/s, ±10°,
                                 as plot_fov.py); "positions" tab: seen from above the ecliptic
                                 with the Parker spirals, and r, lon, lat, source longitudes
                                 (JPL Horizons, cache/horizons_positions.json)
    download original data       the full-resolution AIA frame (JSOC, 1-3 min; then shown at
                                 2048 px instead of the synoptic image, also next time) and the
                                 HRIEUV / HRI-Lya / PHI-HRT files of the step (the file nearest
                                 the time in each shown block), in cache/aia_full/, cache/solo/
    show Solar Orbiter originals those files below the images, each in its own view at full
                                 resolution (HRIEUV, HRI-Lya: asinh stretch; PHI-HRT: B_LOS
                                 ±100 G), zoom / pan, cursor position below; downloaded first if needed

A step is [time - cadence/2, time + cadence/2]. Steps before and after are prepared in the
background. Images are reduced once to a 1024-px north-up PNG with its extent in arcsec
(cache/web/); FOV outlines are computed in closed form (no WCS plotting). Settings and the
last date are kept in browse_web_settings.json.

Run:  ./browse_web.py                     # opens http://localhost:5006/browse
      ./browse_web.py 2025-10-20T12:00 --port 5007 --no-show
"""
import argparse
import html
import json
import threading
import traceback
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import astropy.units as u
import matplotlib
import numpy as np
import pandas as pd
import panel as pn
import sunpy.map
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval
from bokeh.events import Tap
from astropy.time import Time
from bokeh.events import MouseLeave, MouseMove
from bokeh.layouts import column
from bokeh.models import (BoxAnnotation, ColumnDataSource, CustomJS, Div, FactorRange, HoverTool, Label, Range1d,
                          Span)
from bokeh.plotting import figure
from matplotlib.colors import to_hex
from PIL import Image
from sunpy.coordinates import get_earth

import browse as br
import dkist_meta as dm
import plot_fov as pf
import soar_meta as sm
import soops
import spice_meta as sp

OUT_DIR = Path(__file__).resolve().parent
SETTINGS = OUT_DIR / "browse_web_settings.json"
LIGHT = OUT_DIR / "cache" / "web"       # reduced images: PNG + JSON (extent, observer)
LIGHT_N = 1024
FULL_N = 2048                           # full-resolution AIA and the Solar Orbiter originals
BLOS_RANGE = 100                        # G, grey scale of PHI-HRT B_LOS
SC_NAMES = ["Solar Orbiter", "PSP"]     # "positions": these, with their source regions
DEFAULTS = {"wavelength": 171, "cadence": "1h", "show": br.SHOW, "date": "2025-10-20T12:00:00", "positions": False}
ROWS = ["HRIEUV", "HRI-Lya", "PHI-HRT", "SPICE", "FSI", "DKIST"]     # timeline, top down
FSI_COLOR = "#9a9a9a"
FSI_RUN_GAP = pd.Timedelta("1h")
PANEL_PX = 500
BG = "#111111"
_locks, _locks_lock = {}, threading.Lock()


def lock_for(key):
    """One lock per file, so two threads never fetch or reduce the same file at once."""
    with _locks_lock:
        return _locks.setdefault(str(key), threading.Lock())


def load_settings():
    try:
        return DEFAULTS | json.loads(SETTINGS.read_text())
    except (FileNotFoundError, ValueError):
        return dict(DEFAULTS)


def save_settings(s):
    SETTINGS.write_text(json.dumps(s, indent=1))


# ---------------------------------------------------------------- reduced images

def light_image(path, size=LIGHT_N, style="auto"):
    """A FITS image as an RGBA array (row 0 at the bottom), north up, at most size px a
    side, with its extent in arcsec and its observer: dict(rgba, x0, y0, dw, dh, date,
    obs_lat, obs_lon, obs_r, rsun, source). style "auto": the map's colour map, asinh
    stretch; "blos": grey, +-BLOS_RANGE G (PHI-HRT). Made once, kept in cache/web/."""
    path = Path(path)
    tag = f"{size}" if style == "auto" else f"{size}_{style}"
    png, meta = LIGHT / f"{path.name}_{tag}.png", LIGHT / f"{path.name}_{tag}.json"
    with lock_for(png):
        if not (png.is_file() and meta.is_file()):
            m = sunpy.map.Map(path)
            data = np.asarray(m.data, dtype=np.float32)
            if float(m.meta.get("lvl_num", 0)) == 1.0:   # full-resolution AIA level 1: DN -> DN/s
                data = data / np.float32(m.exposure_time.to_value(u.s))
            m = sunpy.map.Map(np.where(np.isfinite(data), data, 0), m.meta)
            if max(m.data.shape) > size:
                f = max(m.data.shape) / size
                m = m.resample([round(m.data.shape[1] / f), round(m.data.shape[0] / f)] * u.pix)
            if not np.allclose(m.rotation_matrix, np.eye(2), atol=1e-4):
                m = m.rotate(order=1, missing=0.0)
            d = m.data
            if style == "blos":
                v = np.clip((d + BLOS_RANGE) / (2 * BLOS_RANGE), 0, 1)
                rgba = (matplotlib.colormaps["gray"](v) * 255).astype(np.uint8)
            else:
                norm = ImageNormalize(d, interval=PercentileInterval(99.7), stretch=AsinhStretch(0.05))
                rgba = (m.cmap(np.clip(norm(d).filled(0), 0, 1)) * 255).astype(np.uint8)
            ny, nx = d.shape
            c = m.pixel_to_world([-0.5, nx - 0.5] * u.pix, [-0.5, ny - 0.5] * u.pix)
            tx, ty = c.Tx.to_value(u.arcsec), c.Ty.to_value(u.arcsec)
            o = m.observer_coordinate
            info = dict(x0=float(tx[0]), y0=float(ty[0]), dw=float(tx[1] - tx[0]), dh=float(ty[1] - ty[0]),
                        date=m.date.isot, obs_lat=float(o.lat.to_value(u.deg)), obs_lon=float(o.lon.to_value(u.deg)),
                        obs_r=float(o.radius.to_value(u.au)), rsun=float(m.rsun_obs.to_value(u.arcsec)),
                        source=path.name)
            LIGHT.mkdir(parents=True, exist_ok=True)
            Image.fromarray(np.flipud(rgba)).save(png, optimize=False)
            meta.write_text(json.dumps(info, indent=1))
        info = json.loads(meta.read_text())
        info["rgba"] = np.ascontiguousarray(np.flipud(np.asarray(Image.open(png).convert("RGBA"))))
    return info


# ---------------------------------------------------------------- geometry (closed form)

def obs(img):
    """(lat, lon, r_au) of the observer of a reduced image."""
    return img["obs_lat"], img["obs_lon"], img["obs_r"]


def seen(lon, lat, img):
    """True where surface points face the observer of img."""
    return pf.cos_seen(lon, lat, obs(img))


KINDS = ("on", "hidden", "off")   # solid, dotted, dashed (plot_fov.project)


def split(x, y, kind):
    """A polyline as one NaN-broken copy per kind (only that kind's runs drawn); a run shares
    its last point with the next, so the line has no gaps."""
    kind = np.asarray(kind)
    nan = np.full(len(x), np.nan)
    out = {}
    for k in KINDS:
        m = kind == k
        m = m | np.r_[False, m[:-1]]
        out[k] = (np.where(m, x, nan), np.where(m, y, nan))
    return out


def surface_on(img, lon, lat):
    """A line on the surface (Stonyhurst lon, lat) drawn on img: x, y and kind (on / hidden)."""
    ok = np.isfinite(lon) & np.isfinite(lat)
    lon, lat = lon[ok], lat[ok]
    x, y, _ = pf.hgs_to_hpc(lon, lat, *obs(img))
    return x, y, np.where(seen(lon, lat, img), "on", "hidden")


def observed_on(tx, ty, src):
    """An outline observed from the image's own vantage point src: kind on (disk) / off."""
    lon, _ = sm.hpc_to_hgs(tx, ty, *src)
    return tx, ty, np.where(np.isfinite(lon), "on", "off")


def limb_of(lat_o, lon_o, r_au, n=181):
    """Stonyhurst lon, lat of the limb of an observer (where its line of sight grazes the surface)."""
    rho = np.arccos(sm.RSUN_KM / (r_au * sm.AU_KM))
    la, lo = np.radians(lat_o), np.radians(lon_o)
    o = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    e1 = np.cross([0, 0, 1], o)
    e1 = e1 / np.linalg.norm(e1) if np.linalg.norm(e1) > 1e-9 else np.array([1.0, 0, 0])
    e2 = np.cross(o, e1)
    phi = np.linspace(0, 2 * np.pi, n)
    p = np.cos(rho) * o[:, None] + np.sin(rho) * (np.cos(phi) * e1[:, None] + np.sin(phi) * e2[:, None])
    return np.degrees(np.arctan2(p[1], p[0])), np.degrees(np.arcsin(np.clip(p[2], -1, 1)))


def grid(img, step=15):
    """Visible lat / lon grid lines of img (NaN-broken), as xs, ys lists."""
    xs, ys = [], []
    lines = [(np.linspace(-180, 180, 361), np.full(361, lat)) for lat in range(-90 + step, 90, step)]
    lines += [(np.full(181, lon), np.linspace(-90, 90, 181)) for lon in range(-180, 180, step)]
    for lon, lat in lines:
        x, y, kind = surface_on(img, lon, lat)
        x, y = split(x, y, kind)["on"]
        xs.append(x)
        ys.append(y)
    return xs, ys


def sub_point(img, lat_o, lon_o):
    """Where an observer's sub-point appears on img: (x, y, front); at the limb if hidden."""
    x, y, _ = pf.hgs_to_hpc(lon_o, lat_o, *obs(img))
    front = bool(seen(lon_o, lat_o, img))
    if not front:
        a = np.arctan2(y, x)
        x, y = img["rsun"] * np.cos(a), img["rsun"] * np.sin(a)
    return float(x), float(y), front


# ---------------------------------------------------------------- one step

class WebStep:
    """Everything one step shows, prepared in a worker thread: the metadata as in
    browse.Step (so browse.summary_lines works on it) plus reduced images and outlines."""

    def __init__(self, t, cadence, wl, store, status, positions=False):
        self.t, self.cadence, self.wl = t, cadence, wl
        half = max(pd.Timedelta(br.CADENCES[cadence]) / 2, pd.Timedelta("1s"))
        self.t0, self.t1 = t - half, t + half
        self.blocks, self.files, fsi = store.window(self.t0, self.t1, status)
        self.fsi_files = fsi
        path = br.full_res_path(t, wl)
        self.full = path is not None
        if path is None:
            status(f"AIA {wl} synoptic for {t:%Y-%m-%d %H:%M} ...")
            with lock_for(f"aia{wl}{t.round('2min')}"):
                path = pf.aia_synoptic(t, wl)
        self.aia = light_image(path, FULL_N if self.full else LIGHT_N)
        try:
            self.dkist = dm.datasets(self.t0, self.t1)
            if len(self.dkist):
                self.dkist = self.dkist.dropna(subset=["x0"])
        except Exception as e:   # network trouble: no DKIST, the rest is shown
            status(f"DKIST: {e}")
            self.dkist = pd.DataFrame()
        self.fsi = None
        if len(self.blocks) and len(fsi):
            g = fsi[fsi["product"] == "fsi174"]
            if len(g):
                dt = (g.start - t).abs()
                k = dt.idxmin()
                if dt[k] <= br.FSI_WINDOW:
                    status(f"FSI 174 {g.start[k]:%Y-%m-%d %H:%M} ...")
                    with lock_for(g.filename[k]):
                        p = sm.download(g.filename[k], pf.FSI_CACHE)
                    self.fsi = light_image(p)
        # campaigns (ESA SOOP pages) and SPICE studies, as browse.Step
        self.soop_pages, self.soop_runs, self.spice = {}, {}, {}
        for soop in dict.fromkeys(self.blocks.soop) if len(self.blocks) else []:
            if str(soop).lower() == "none":
                continue
            status(f"SOOP {soop} ...")
            self.soop_pages[soop] = soops.page(soop)
            r = soops.runs(soop, self.t0, self.t1, margin="0D")
            self.soop_runs[soop] = r.iloc[0] if len(r) else None
        if len(self.blocks) and "study" in self.files:
            spf = self.files[self.files.block.isin(self.blocks.block[self.blocks["product"] == "spice"])]
            for study, g in spf.dropna(subset=["study"]).groupby("study"):
                in_step = g[(g.start >= self.t0) & (g.start <= self.t1)]
                if len(in_step) or (g.start - t).abs().min() <= pd.Timedelta(br.CADENCES[cadence]):
                    status(f"SPICE study {study} ...")
                    self.spice[study] = sp.study(study, g.filename.iloc[0])
        self.outlines = self.make_outlines()
        self.earth = get_earth(Time(t.to_pydatetime()))
        self.coords = {}
        if positions:   # JPL Horizons, cached in cache/horizons_positions.json
            status("spacecraft positions ...")
            self.coords = {n: c for n, c in pf.positions(t, SC_NAMES, {}).items() if c is not None}

    def make_outlines(self):
        """Outlines on AIA and FSI: list of dict(show, label, color, x, y, kind, tip); kind per
        point: on (solid), hidden (dotted, behind the disk), off (dashed, off the limb)."""
        out = {"aia": [], "fsi": []}
        for _, b in (self.blocks.iterrows() if len(self.blocks) else []):
            if not np.isfinite(b.get("c0x", np.nan)):
                continue
            tx, ty = pf.box_outline(b[sm.CORNER_COLS].to_numpy(float).reshape(4, 2))
            src = (b.solo_lat, b.solo_lon, b.solo_r)
            label = pf.LABELS[b["product"]]
            tip = dict(name=f"{label} ({b['product']})",
                       when=f"{b.start:%Y-%m-%d %H:%M} - {b.end:%m-%d %H:%M} ({br.fmt_td(b.end - b.start)})",
                       cadence=(f"one every {br.fmt_s(b.cadence_s)}, " if int(b.n) > 1 else "") + f"{int(b.n)} files",
                       where=br.fov_text(b), soop=str(b.soop))
            base = dict(show=label, label=label, color=pf.COLORS[b["product"]], tip=tip)
            x, y, kind = pf.project(tx, ty, src, obs(self.aia), coronal=b["product"] in pf.CORONAL)
            out["aia"].append(base | dict(x=x, y=y, kind=kind))
            if self.fsi is not None:
                x, y, kind = observed_on(tx, ty, src)
                out["fsi"].append(base | dict(x=x, y=y, kind=kind))
        dk = self.dkist
        if len(dk):
            key = dk[["x0", "y0", "x1", "y1"]].round(-1).astype(str).agg(",".join, axis=1)
            dk = dk.assign(key=dk.instrument + key).drop_duplicates("key")
            earth = obs(self.aia)
            for _, r in dk.iterrows():
                inst = str(r.instrument).upper()
                corners = np.array([[r.x0, r.y0], [r.x1, r.y0], [r.x1, r.y1], [r.x0, r.y1]])
                tx, ty = pf.box_outline(corners)
                tip = dict(name=f"DKIST {br.DKIST_NAME.get(inst, inst)}",
                           when=f"{r.start:%Y-%m-%d %H:%M} - {r.end:%m-%d %H:%M}",
                           cadence=f"{r.get('frameCount', '?')} frames", where=f"lines: {r.lines}",
                           soop=f"proposal {r.proposal}")
                base = dict(show="DKIST", label=f"DKIST {br.DKIST_NAME.get(inst, inst)}",
                            color=pf.DKIST_COLORS.get(inst, pf.DKIST_OTHER), tip=tip)
                x, y, kind = observed_on(tx, ty, earth)
                out["aia"].append(base | dict(x=x, y=y, kind=kind))
                if self.fsi is not None:   # off-limb DKIST (corona) goes up off the limb on FSI too
                    x, y, kind = pf.project(tx, ty, earth, obs(self.fsi), coronal=True)
                    out["fsi"].append(base | dict(x=x, y=y, kind=kind))
        return out


# ---------------------------------------------------------------- timeline data

def timeline_month(month):
    """Bars of one month (as loaded by browse.MetaStore): DataFrame row, start, end, color,
    name, info."""
    rows = []
    b = month["blocks"]
    for _, r in (b.iterrows() if len(b) else []):
        label = pf.LABELS.get(r["product"])
        if label is None:
            continue
        cad = f"one every {br.fmt_s(r.cadence_s)}" if int(r.n) > 1 else "single file"
        where = f"centre lat {r.lat:+.0f}°" if np.isfinite(r.get("lat", np.nan)) else "centre off disk / unknown"
        rows.append(dict(row=label, start=r.start, end=r.end, color=pf.COLORS[r["product"]],
                         name=f"{label} ({r['product']})", info=f"{int(r.n)} files, {cad}; {where}; SOOP {r.soop}"))
    f = month["fsi"]
    if len(f):
        g = f[f["product"] == "fsi174"].sort_values("start")
        if len(g):
            run = (g.start.diff() > FSI_RUN_GAP).cumsum()
            for _, x in g.groupby(run):
                cad = x.start.diff().dt.total_seconds().median()
                rows.append(dict(row="FSI", start=x.start.iloc[0], end=x.start.iloc[-1], color=FSI_COLOR,
                                 name="FSI 174", info=f"{len(x)} images, one every {br.fmt_s(cad)}"))
    return rows


def timeline_dkist(m):
    end = m + pd.offsets.MonthBegin(1)
    try:
        d = dm.datasets(m, end)
    except Exception:
        return []
    rows = []
    for _, r in (d.iterrows() if len(d) else []):
        inst = str(r.instrument).upper()
        rows.append(dict(row="DKIST", start=r.start, end=r.end, color=pf.DKIST_COLORS.get(inst, pf.DKIST_OTHER),
                         name=f"DKIST {br.DKIST_NAME.get(inst, inst)}",
                         info=f"proposal {r.proposal}; lines {r.lines}"))
    return rows


# ---------------------------------------------------------------- spacecraft on the images, positions plot

def observer_layers(img, bodies, connect):
    """What to draw on img for observers and source regions:
    bodies   [(name, lat, lon, r_au, color)]: sub-point (x; o at the limb when hidden) and limb
    connect  [(name, SkyCoord, color)]: connectivity region (Parker spiral to 2.5 Rsun,
             plot_fov.V_RANGE km/s, +-plot_fov.MARGIN deg) and the speed line along its latitude
    -> data for the sources mark, limb, sc_on, sc_hidden."""
    mark = dict(x=[], y=[], color=[], text=[], marker=[])
    limb = dict(xs=[], ys=[], color=[])
    lines = {"on": dict(xs=[], ys=[], color=[]), "hidden": dict(xs=[], ys=[], color=[])}
    for name, lat, lon, r, color in bodies:
        x, y, front = sub_point(img, lat, lon)
        for k, v in dict(x=x, y=y, color=color, marker="x" if front else "circle",
                         text=name if front else f"{name} (far side)").items():
            mark[k].append(v)
        lx, ly, kind = surface_on(img, *limb_of(lat, lon, r))
        sx, sy = split(lx, ly, kind)["on"]
        limb["xs"].append(sx)
        limb["ys"].append(sy)
        limb["color"].append(color)
    for name, c, color in connect:
        fast, slow, lon, lat = pf.connectivity(c, pf.V_RANGE, pf.MARGIN)
        line_lon = np.linspace(fast, slow, 40)
        for lo, la in ((lon, lat), (line_lon, np.full_like(line_lon, c.lat.to_value(u.deg)))):
            x, y, kind = surface_on(img, lo, la)
            parts = split(x, y, kind)
            for k in lines:
                lines[k]["xs"].append(parts[k][0])
                lines[k]["ys"].append(parts[k][1])
                lines[k]["color"].append(color)
        x, y, kind = surface_on(img, lon, lat)
        top = 2 * len(lon) // 4 + len(lon) // 8   # middle of the north edge
        if len(x) > top and kind[top] == "on":
            for k, v in dict(x=x[top], y=y[top], color=color, marker="dot", text=f"{name} source").items():
                mark[k].append(v)
    return dict(mark=mark, limb=limb, sc_on=lines["on"], sc_hidden=lines["hidden"])


def positions_figure():
    """Spacecraft seen from above the ecliptic (Stonyhurst, au; Earth at lon 0, to the right)."""
    p = figure(frame_width=430, frame_height=430, match_aspect=True, tools="pan,wheel_zoom,reset",
               active_scroll="wheel_zoom", background_fill_color=BG, border_fill_color=BG, outline_line_color=None,
               x_range=Range1d(-1.1, 1.1), y_range=Range1d(-1.1, 1.1), title="positions (Stonyhurst, au; Earth at 0°)")
    p.title.text_color = "#eeeeee"
    p.grid.visible = False
    for ax in (p.xaxis, p.yaxis):
        ax.major_label_text_color = "#cccccc"
    ph = np.linspace(0, 2 * np.pi, 200)
    p.multi_line([list(r * np.cos(ph)) for r in (0.25, 0.5, 0.75, 1.0)],
                 [list(r * np.sin(ph)) for r in (0.25, 0.5, 0.75, 1.0)], line_color="#444444")
    p.text([0.25, 0.5, 0.75, 1.0], [0, 0, 0, 0], text=["0.25", "0.5", "0.75", "1 au"], text_color="#777777",
           text_font_size="9px", y_offset=12)
    p.scatter([0], [0], size=14, color="orange")
    src = dict(band=ColumnDataSource(dict(xs=[], ys=[], color=[])),
               spiral=ColumnDataSource(dict(xs=[], ys=[], color=[], dash=[])),
               ray=ColumnDataSource(dict(xs=[], ys=[], color=[])),
               body=ColumnDataSource(dict(x=[], y=[], color=[], name=[], marker=[], info=[])))
    p.patches("xs", "ys", source=src["band"], fill_color="color", fill_alpha=0.15, line_alpha=0)
    p.multi_line("xs", "ys", source=src["spiral"], line_color="color", line_dash="dashed", line_width=1)
    p.multi_line("xs", "ys", source=src["ray"], line_color="color", line_alpha=0.5)
    r = p.scatter("x", "y", source=src["body"], marker="marker", size=10, color="color", line_color="black")
    p.text("x", "y", text="name", source=src["body"], text_color="color", text_font_size="11px", x_offset=7, y_offset=-4)
    p.add_tools(HoverTool(renderers=[r], tooltips=[("", "@name"), ("", "@info")]))
    return p, src


def positions_data(st):
    """Data of positions_figure() and an HTML table, for a step with coords."""
    from sunpy.coordinates.sun import L0 as carrington_L0
    bodies = {"Earth": st.earth} | st.coords
    band, spiral, ray, body, rows = ([] for _ in range(5))
    rmax = 0
    car = carrington_L0(Time(st.t.to_pydatetime())).to_value(u.deg)
    for n, c in bodies.items():
        color = pf.EARTH_COLOR if n == "Earth" else pf.SC_COLORS.get(n, "white")
        lon, lat, r = c.lon.wrap_at(180 * u.deg).to_value(u.deg), c.lat.to_value(u.deg), c.radius.to_value(u.au)
        rmax = max(rmax, r)
        x, y = r * np.cos(np.radians(lon)), r * np.sin(np.radians(lon))
        sep = np.degrees(np.arccos(np.clip(sm.cos_angle(lon, lat, 0.0, st.earth.lat.to_value(u.deg)), -1, 1)))
        info = f"r {r:.3f} au ({c.radius.to_value(u.R_sun):.0f} Rsun), lon {lon:+.1f}°, lat {lat:+.1f}°, {sep:.0f}° from Earth"
        if n != "Earth":
            rr = np.linspace(pf.R_SS.to_value(u.au), r, 200) * u.au
            th = [np.radians(lon + pf.spiral_shift(c.radius, rr, v)) for v in pf.V_RANGE]
            xs = [rr.value * np.cos(t_) for t_ in th]
            ys = [rr.value * np.sin(t_) for t_ in th]
            band.append((np.r_[xs[0], xs[1][::-1]], np.r_[ys[0], ys[1][::-1]], color))
            spiral += [(xs[k], ys[k], color, d) for k, d in ((0, "dashed"), (1, "dotted"))]
            fast, slow, _, _ = pf.connectivity(c, pf.V_RANGE, pf.MARGIN)
            info += (f"; connects lon {fast:+.0f}..{slow:+.0f} (Carrington {(fast + car) % 360:.0f}.."
                     f"{(slow + car) % 360:.0f}) ±{pf.MARGIN}° for {pf.V_RANGE[0]}-{pf.V_RANGE[1]} km/s")
        ray.append(([0, x], [0, y], color))
        body.append((x, y, color, n, "circle" if n == "Earth" else "diamond", info))
        rows.append(f'<div style="color:{to_hex(color)}"><b>{n}</b>: {info}</div>')
    data = dict(band=dict(xs=[b[0] for b in band], ys=[b[1] for b in band], color=[b[2] for b in band]),
                spiral=dict(xs=[q[0] for q in spiral], ys=[q[1] for q in spiral], color=[q[2] for q in spiral],
                            dash=[q[3] for q in spiral]),
                ray=dict(xs=[q[0] for q in ray], ys=[q[1] for q in ray], color=[q[2] for q in ray]),
                body={k: [b[i] for b in body] for i, k in enumerate(("x", "y", "color", "name", "marker", "info"))})
    note = (f"source regions: ballistic Parker spiral to 2.5 Rsun, {pf.V_RANGE[0]}-{pf.V_RANGE[1]} km/s "
            f"(dashed fast, dotted slow), ±{pf.MARGIN}°, no PFSS: in a coronal hole the photospheric footpoint can lie elsewhere")
    table = ('<div style="font-size:12px;color:#ccc">' + "".join(rows) +
             f'<div style="color:#888;margin-top:4px">{note}</div></div>')
    return data, 1.15 * rmax, table


# ---------------------------------------------------------------- original data

SOLO_FILES = br.SOLO_FILES
ORIGINAL_PRODUCTS = ("hri174", "hriopn", "hrilya", "phi-hrt")   # SPICE rasters are large: not fetched


def solo_rows(st, products):
    """For each shown HRIEUV / HRI-Lya / PHI-HRT block of the step, its file nearest the time."""
    out = []
    b = st.blocks
    if not len(b):
        return out
    for _, blk in b[b["product"].isin([p for p in products if p in ORIGINAL_PRODUCTS])].iterrows():
        g = st.files[st.files.block == blk.block]
        if len(g):
            out.append(g.iloc[(g.start - st.t).abs().argmin()])
    return out


def prepare_originals(rows, status):
    """Download (once, cache/solo/<product>/) and reduce the Solar Orbiter files:
    [(product, filename, image)]."""
    out = []
    for r in rows:
        status(f"Solar Orbiter original: {r.filename} ...")
        with lock_for(r.filename):
            path = sm.download(r.filename, SOLO_FILES / r["product"])
        style = "blos" if r["product"] == "phi-hrt" else "auto"
        out.append((r["product"], r.filename, light_image(path, FULL_N, style)))
    return out


def original_figure(product, filename, img):
    """One Solar Orbiter file, in its own view (helioprojective arcsec, north up)."""
    unit = f"B_LOS ±{BLOS_RANGE} G" if product == "phi-hrt" else "asinh stretch"
    p = figure(title=f"{pf.LABELS[product]} {img['date'][:19].replace('T', ' ')}  ({unit})",
               frame_width=560, frame_height=560, match_aspect=True,
               x_range=Range1d(img["x0"], img["x0"] + img["dw"]), y_range=Range1d(img["y0"], img["y0"] + img["dh"]),
               tools="pan,wheel_zoom,box_zoom,reset,save", active_scroll="wheel_zoom",
               background_fill_color="black", border_fill_color=BG, outline_line_color=None)
    p.title.text_color = pf.COLORS[product]
    for ax in (p.xaxis, p.yaxis):
        ax.major_label_text_color = "#cccccc"
    ny, nx = img["rgba"].shape[:2]
    src = ColumnDataSource(dict(image=[img["rgba"].view(np.uint32).reshape(ny, nx)], x=[img["x0"]], y=[img["y0"]],
                                dw=[img["dw"]], dh=[img["dh"]], **image_meta(img)))
    p.image_rgba(image="image", x="x", y="y", dw="dw", dh="dh", source=src)
    p.grid.visible = False
    p.add_layout(Label(x=5, y=5, x_units="screen", y_units="screen", text=filename, text_color="#888888",
                       text_font_size="10px"))
    return with_readout(p, src, 640)


# ---------------------------------------------------------------- the page

class Store(br.MetaStore):
    """browse.MetaStore, but the months next to a step only for FSI (a step near the end of a
    month needs the next month's FSI images, not all its Solar Orbiter metadata)."""

    def __init__(self):
        super().__init__()
        self.fsi_months = {}

    def fsi_month(self, m):
        with self.lock:
            if m not in self.fsi_months:
                end = m + pd.offsets.MonthBegin(1) - pd.Timedelta("1us")
                self.fsi_months[m] = sm.files(m, end, ["fsi174", "fsi304"], headers=False)
            return self.fsi_months[m]

    def window(self, t0, t1, status=print):
        parts = [self.month(m, status) for m in sm.months(t0, t1)]
        nonempty = lambda key: [p[key] for p in parts if len(p[key])]
        b = pd.concat(nonempty("blocks"), ignore_index=True) if nonempty("blocks") else pd.DataFrame()
        f = pd.concat(nonempty("files"), ignore_index=True) if nonempty("files") else pd.DataFrame()
        fs = [self.months[m]["fsi"] if m in self.months else self.fsi_month(m)
              for m in sm.months(t0 - br.FSI_WINDOW, t1 + br.FSI_WINDOW)]
        fs = [x for x in fs if len(x)]
        fsi = pd.concat(fs, ignore_index=True) if fs else pd.DataFrame()
        if len(b):
            b = b[(b.start <= t1) & (b.end >= t0)]
        return b, f, fsi


STORE = Store()                 # shared by all sessions
POOL = ThreadPoolExecutor(max_workers=4)


# The cursor position under an image, written in the Div below it, in the browser (no server round trip):
# helioprojective x, y and Stonyhurst lon, lat and mu (closed form as soar_meta.hpc_to_hgs), or off the
# disk the distance from the Sun's centre in Rsun; the observer is in the image's source (obs_*, rsun)
READOUT_JS = """
const d = src.data;
if (!d.obs_r || d.obs_r.length == 0) { div.text = ""; return; }
const X = cb_obj.x, Y = cb_obj.y;
const rad = Math.PI / 180, D = d.obs_r[0] * %(au)r, R = %(rsun)r;
const B0 = d.obs_lat[0] * rad, L0 = d.obs_lon[0] * rad;
const tx = X / 3600 * rad, ty = Y / 3600 * rad;
const q = D * Math.cos(ty) * Math.cos(tx), disc = q * q - D * D + R * R;
const hpc = "<b>Helioprojective</b> x " + X.toFixed(1) + "\u2033, y " + Y.toFixed(1) + "\u2033";
let hgs;
if (disc < 0) {
    hgs = "off the disk, " + (Math.hypot(X, Y) / d.rsun[0]).toFixed(3) + " R\u2609 from centre (plane of sky)";
} else {
    const s = q - Math.sqrt(disc);
    const x = s * Math.cos(ty) * Math.sin(tx), y = s * Math.sin(ty), z = D - s * Math.cos(ty) * Math.cos(tx);
    const lat = Math.asin(Math.max(-1, Math.min(1, (y * Math.cos(B0) + z * Math.sin(B0)) / R))) / rad;
    let lon = (L0 + Math.atan2(x, z * Math.cos(B0) - y * Math.sin(B0))) / rad;
    lon = ((lon + 540) %% 360) - 180;
    // mu: surface normal . unit vector from the point to the observer (exact for a finite distance)
    const mu = (-x * Math.sin(tx) * Math.cos(ty) - y * Math.sin(ty) + z * Math.cos(ty) * Math.cos(tx)) / R;
    hgs = "lon " + lon.toFixed(2) + "\u00b0, lat " + lat.toFixed(2) + "\u00b0, \u03bc " + mu.toFixed(3);
}
div.text = hpc + " &nbsp;&nbsp; <b>Stonyhurst</b> " + hgs;
""" % dict(au=sm.AU_KM, rsun=sm.RSUN_KM)


def with_readout(p, src, width):
    """The figure with a line below it giving the cursor position (see READOUT_JS)."""
    div = Div(text="", width=width, height=22,
              styles={"color": "#dddddd", "font-family": "monospace", "font-size": "12px", "padding-left": "8px"})
    p.js_on_event(MouseMove, CustomJS(args=dict(src=src, div=div), code=READOUT_JS))
    p.js_on_event(MouseLeave, CustomJS(args=dict(div=div), code="div.text = '';"))
    return column(p, div)


def image_meta(img):
    """Image source columns: the observer, for the cursor readout (READOUT_JS)."""
    return dict(obs_lat=[img["obs_lat"]], obs_lon=[img["obs_lon"]], obs_r=[img["obs_r"]], rsun=[img["rsun"]])


def image_figure(title):
    p = figure(title=title, frame_width=PANEL_PX, frame_height=PANEL_PX, match_aspect=True,
               x_range=Range1d(-1200, 1200), y_range=Range1d(-1200, 1200),
               tools="pan,wheel_zoom,box_zoom,reset,save", active_scroll="wheel_zoom",
               background_fill_color="black", border_fill_color=BG, outline_line_color=None)
    p.xaxis.axis_label = "x (arcsec)"
    p.yaxis.axis_label = "y (arcsec)"
    for ax in (p.xaxis, p.yaxis):
        ax.axis_label_text_color = ax.major_label_text_color = "#cccccc"
    p.title.text_color = "#eeeeee"
    p.grid.visible = False
    src = dict(img=ColumnDataSource(dict(image=[], x=[], y=[], dw=[], dh=[], obs_lat=[], obs_lon=[], obs_r=[], rsun=[])),
               grid=ColumnDataSource(dict(xs=[], ys=[])),
               on=ColumnDataSource(dict(xs=[], ys=[], color=[])),
               hidden=ColumnDataSource(dict(xs=[], ys=[], color=[])),
               off=ColumnDataSource(dict(xs=[], ys=[], color=[])),
               patch=ColumnDataSource(dict(xs=[], ys=[], color=[], label=[], name=[], when=[], cadence=[],
                                           where=[], soop=[])),
               limb=ColumnDataSource(dict(xs=[], ys=[], color=[])),
               sc_on=ColumnDataSource(dict(xs=[], ys=[], color=[])),
               sc_hidden=ColumnDataSource(dict(xs=[], ys=[], color=[])),
               mark=ColumnDataSource(dict(x=[], y=[], color=[], text=[], marker=[])))
    p.image_rgba(image="image", x="x", y="y", dw="dw", dh="dh", source=src["img"])
    p.multi_line("xs", "ys", source=src["grid"], line_color="white", line_alpha=0.22, line_width=0.6)
    patch = p.patches("xs", "ys", source=src["patch"], fill_color="color", fill_alpha=0.07, line_alpha=0,
                      legend_field="label")
    p.multi_line("xs", "ys", source=src["limb"], line_color="color", line_width=1.2, line_dash="dashed")
    p.multi_line("xs", "ys", source=src["sc_on"], line_color="color", line_width=2)
    p.multi_line("xs", "ys", source=src["sc_hidden"], line_color="color", line_width=2, line_dash="dotted")
    for key, dash in (("on", "solid"), ("hidden", "dotted"), ("off", "dashed")):
        p.multi_line("xs", "ys", source=src[key], line_color="black", line_width=4, line_alpha=0.6, line_dash=dash)
        r = p.multi_line("xs", "ys", source=src[key], line_color="color", line_width=2, line_dash=dash)
        src[key + "_r"] = r
    p.scatter("x", "y", source=src["mark"], marker="marker", size=12, line_color="color", fill_alpha=0, line_width=2.5)
    p.text("x", "y", text="text", source=src["mark"], text_color="color", text_font_size="11px",
           x_offset=8, y_offset=-4)
    p.add_tools(HoverTool(renderers=[patch], tooltips=[("", "@name"), ("time", "@when"), ("files", "@cadence"),
                                                       ("FOV", "@where"), ("SOOP", "@soop")]))
    p.legend.location = "bottom_left"
    p.legend.background_fill_alpha = 0.5
    p.legend.background_fill_color = "black"
    p.legend.label_text_color = "white"
    p.legend.label_text_font_size = "10px"
    note = Label(x=10, y=PANEL_PX // 2, x_units="screen", y_units="screen", text="", text_color="#888888",
                 text_font_size="13px")
    p.add_layout(note)
    src["note"] = note
    src["layout"] = with_readout(p, src["img"], PANEL_PX + 80)
    return p, src


class App:
    def __init__(self, t=None):
        self.s = load_settings()
        self.t = br.parse_time(t or self.s["date"])
        self.steps = {}           # key -> future
        self.pending = None
        self.status_text = ""
        self.months_drawn = set()
        self.tl_rows = {}         # month -> bars
        self.first = {"aia": True, "fsi": True}
        self.orig_job = self.solo_job = None
        self.slow = ThreadPoolExecutor(max_workers=1)    # original data downloads (JSOC: minutes)

        # controls
        self.date = pn.widgets.TextInput(name="date (UTC, Enter)", value=f"{self.t:%Y-%m-%d %H:%M:%S}", width=190)
        self.back = pn.widgets.Button(name="◀ back", width=80)
        self.fwd = pn.widgets.Button(name="forward ▶", width=90)
        self.cad = pn.widgets.Select(name="cadence", options=list(br.CADENCES), value=self.s["cadence"], width=90)
        self.wl = pn.widgets.Select(name="AIA (Å)", options=[str(w) for w in br.WAVELENGTHS],
                                    value=str(self.s["wavelength"]), width=90)
        self.show = pn.widgets.CheckButtonGroup(name="show", options=br.SHOW, value=list(self.s["show"]))
        self.orig_btn = pn.widgets.Button(name="download original data", button_type="primary", width=180,
                                          description="full-resolution AIA (JSOC, 1-3 min) and the Solar Orbiter "
                                                      "files nearest the time; kept in cache/")
        self.pos_btn = pn.widgets.Toggle(name="positions", value=bool(self.s.get("positions")), width=100)
        self.solo_btn = pn.widgets.Button(name="show Solar Orbiter originals", width=210,
                                          description="the HRIEUV / HRI-Lya / PHI-HRT files nearest the time "
                                                      "(downloaded first if needed)")
        self.orig_btn.on_click(lambda _: self.download_originals())
        self.solo_btn.on_click(lambda _: self.show_originals())
        self.pos_btn.param.watch(self.positions_toggled, "value")
        self.status = pn.pane.HTML("", styles={"color": "#ff8080"}, width=450)
        self.date.param.watch(self.date_typed, "value")
        self.back.on_click(lambda _: self.go(self.t - self.cadence()))
        self.fwd.on_click(lambda _: self.go(self.t + self.cadence()))
        self.cad.param.watch(self.setting_changed, "value")
        self.wl.param.watch(self.setting_changed, "value")
        self.show.param.watch(lambda _: self.show_changed(), "value")

        # timeline
        self.tl_src = ColumnDataSource(dict(row=[], start=[], end=[], color=[], name=[], info=[]))
        tl = figure(x_axis_type="datetime", y_range=FactorRange(*reversed(ROWS)), height=210,
                    sizing_mode="stretch_width", tools="xpan,xwheel_zoom,reset", active_scroll="xwheel_zoom",
                    active_drag="xpan", background_fill_color=BG, border_fill_color=BG, outline_line_color=None,
                    x_range=Range1d(self.t - pd.Timedelta("3D"), self.t + pd.Timedelta("3D")))
        bars = tl.hbar(y="row", left="start", right="end", height=0.6, fill_color="color", line_color="color",
                       line_width=2, source=self.tl_src)
        tl.add_tools(HoverTool(renderers=[bars], tooltips=[("", "@name"), ("from", "@start{%F %H:%M}"),
                                                           ("to", "@end{%F %H:%M}"), ("", "@info")],
                               formatters={"@start": "datetime", "@end": "datetime"}))
        self.band = BoxAnnotation(fill_color="white", fill_alpha=0.15, line_color="white", line_alpha=0.6)
        self.now = Span(dimension="height", line_color="white", line_width=1)
        tl.add_layout(self.band)
        tl.add_layout(self.now)
        for ax in (tl.xaxis, tl.yaxis):
            ax.major_label_text_color = "#cccccc"
        tl.grid.grid_line_alpha = 0.15
        tl.on_event(Tap, self.timeline_tap)
        self.tl = tl

        # images and tabs
        self.fig_aia, self.src_aia = image_figure("SDO/AIA")
        self.fig_fsi, self.src_fsi = image_figure("SolO/EUI FSI 174 Å")
        self.summary = pn.pane.HTML("", sizing_mode="stretch_both")
        self.files_tab = pn.widgets.Tabulator(pd.DataFrame(), disabled=True, show_index=False, sizing_mode="stretch_both")
        self.dkist_tab = pn.widgets.Tabulator(pd.DataFrame(), disabled=True, show_index=False, sizing_mode="stretch_both")
        self.title = pn.pane.HTML("", styles={"font-size": "15px", "color": "#eeeeee"})
        self.fig_pos, self.src_pos = positions_figure()
        self.pos_table = pn.pane.HTML("", sizing_mode="stretch_width")
        self.pos_note = pn.pane.HTML('<span style="color:#888">press "positions" to see where the spacecraft are</span>')
        self.tabs = pn.Tabs(("summary", self.summary), ("Solar Orbiter files", self.files_tab), ("DKIST", self.dkist_tab),
                            ("positions", pn.Column(self.pos_note, pn.pane.Bokeh(self.fig_pos), self.pos_table)),
                            sizing_mode="stretch_both", min_width=480)
        controls = pn.Row(self.date, self.back, self.fwd, self.cad, self.wl,
                          pn.Column(pn.pane.HTML("show", styles={"color": "#aaaaaa", "font-size": "12px"}),
                                    self.show, margin=0),
                          pn.Column(pn.Spacer(height=18), pn.Row(self.pos_btn, self.orig_btn, self.solo_btn), margin=0),
                          sizing_mode="stretch_width")
        self.close_btn = pn.widgets.Button(name="close", width=80)
        self.close_btn.on_click(lambda _: setattr(self.orig_area, "visible", False))
        self.orig_area = pn.Column(sizing_mode="stretch_width", visible=False)
        self.layout = pn.Column(controls, pn.Row(self.title, self.status), pn.pane.Bokeh(tl, sizing_mode="stretch_width"),
                                pn.Row(pn.pane.Bokeh(self.src_aia["layout"]), pn.pane.Bokeh(self.src_fsi["layout"]), self.tabs,
                                       sizing_mode="stretch_width", height=PANEL_PX + 120),
                                self.orig_area, sizing_mode="stretch_width", styles={"background": BG})
        self.go(self.t)

    # ------------------------------------------------ controls
    def cadence(self):
        return pd.Timedelta(br.CADENCES[self.s["cadence"]])

    def set_status(self, text):
        self.status_text = text      # from worker threads; shown by poll()

    def date_typed(self, event):
        if event.new == f"{self.t:%Y-%m-%d %H:%M:%S}":
            return
        try:
            self.go(br.parse_time(event.new))
        except ValueError:
            self.set_status(f"can't read the date '{html.escape(event.new)}'")

    def setting_changed(self, _):
        self.s["cadence"], self.s["wavelength"] = self.cad.value, int(self.wl.value)
        save_settings(self.s)
        self.go(self.t)

    def show_changed(self):
        self.s["show"] = list(self.show.value)
        save_settings(self.s)
        if self.step is not None:
            self.draw()

    def timeline_tap(self, event):
        self.go(pd.Timestamp(event.x, unit="ms").floor("s"))

    def positions_toggled(self, event):
        self.s["positions"] = bool(event.new)
        save_settings(self.s)
        if event.new:
            self.tabs.active = 3
        self.go(self.t)

    def shown_products(self):
        return [p for n in self.s["show"] for p in br.PRODUCT_OF.get(n, [])]

    def download_originals(self):
        """Full-resolution AIA and the Solar Orbiter files of the step, in the background."""
        if self.orig_job is not None:
            self.set_status("original data: already downloading ...")
            return
        st = self.step
        if st is None:
            return
        rows = solo_rows(st, self.shown_products())
        self.set_status(f"original data: AIA {st.wl} full resolution (JSOC, 1-3 min) + {len(rows)} Solar Orbiter file(s) ...")
        self.orig_job = self.slow.submit(br.fetch_originals, st.t, st.wl, st.full, rows, self.set_status)

    def show_originals(self):
        """The Solar Orbiter files nearest the time, below the images (downloaded if needed)."""
        if self.solo_job is not None:
            self.set_status("Solar Orbiter originals: already preparing ...")
            return
        st = self.step
        if st is None:
            return
        rows = solo_rows(st, self.shown_products())
        if not rows:
            self.set_status("no HRIEUV / HRI-Lya / PHI-HRT observation in this step (shown products)")
            return
        self.solo_job = (st.t, POOL.submit(prepare_originals, rows, self.set_status))

    # ------------------------------------------------ stepping
    step = None

    def key(self, t):
        return (t, self.s["cadence"], self.s["wavelength"], bool(self.s.get("positions")))

    def submit(self, t):
        k = self.key(t)
        if k not in self.steps:
            def status(text):   # only the step being waited for reports
                if self.pending is not None and self.pending[0] == k:
                    self.set_status(text)
            self.steps[k] = POOL.submit(WebStep, t, self.s["cadence"], self.s["wavelength"], STORE, status,
                                        bool(self.s.get("positions")))
        return k, self.steps[k]

    def go(self, t):
        self.t = t
        self.s["date"] = f"{t:%Y-%m-%dT%H:%M:%S}"
        save_settings(self.s)
        if self.date.value != f"{t:%Y-%m-%d %H:%M:%S}":
            self.date.value = f"{t:%Y-%m-%d %H:%M:%S}"
        self.pending = self.submit(t)
        self.set_status(f"loading {t:%Y-%m-%d %H:%M:%S} ...")
        half = max(self.cadence() / 2, pd.Timedelta("1s"))
        self.band.left, self.band.right = (t - half).timestamp() * 1000, (t + half).timestamp() * 1000
        self.now.location = t.timestamp() * 1000
        r = self.tl.x_range
        lo, hi = as_time(r.start), as_time(r.end)
        if not (lo <= t <= hi):
            w = hi - lo
            r.start, r.end = t - w / 2, t + w / 2

    def poll(self):
        if self.status.object != self.status_text:
            self.status.object = self.status_text
        if set(STORE.months) - self.months_drawn:
            self.draw_timeline()
        if self.orig_job is not None and self.orig_job.done():
            job, self.orig_job = self.orig_job, None
            try:
                t, paths = job.result()
                self.set_status(f"original data of {t:%Y-%m-%d %H:%M:%S} downloaded ({len(paths)} Solar Orbiter "
                                f"file(s)); \"show Solar Orbiter originals\" to see them")
                self.steps = {k: v for k, v in self.steps.items() if k[0] != t}   # again, now full resolution
                if t == self.t:
                    self.go(t)
            except Exception as e:
                traceback.print_exc()
                self.set_status(f"original data failed: {html.escape(str(e))}")
        if self.solo_job is not None and self.solo_job[1].done():
            (t, job), self.solo_job = self.solo_job, None
            try:
                items = job.result()
            except Exception as e:
                traceback.print_exc()
                self.set_status(f"Solar Orbiter originals failed: {html.escape(str(e))}")
            else:
                self.set_status("")
                head = pn.pane.HTML(f'<b style="color:#eee">Solar Orbiter files nearest to {t:%Y-%m-%d %H:%M:%S}</b>'
                                    '<span style="color:#888"> &nbsp; own view, north up; wheel zoom, drag pan, '
                                    'cursor position below each image</span>')
                self.orig_area.objects = [pn.Row(head, self.close_btn),
                                          pn.Row(*[pn.pane.Bokeh(original_figure(*it)) for it in items])]
                self.orig_area.visible = True
        if self.pending is None or not self.pending[1].done():
            return
        k, fut = self.pending
        self.pending = None
        try:
            self.step = fut.result()
        except Exception as e:
            traceback.print_exc()
            self.set_status(f"failed: {html.escape(str(e))}")
            self.steps.pop(k, None)
            return
        self.set_status("")
        self.draw()
        keep = {self.key(self.t + d * self.cadence()) for d in (-1, 0, 1)}
        for old in [x for x in self.steps if x not in keep]:
            self.steps.pop(old)
        for d in (1, -1):            # prepare the neighbours
            self.submit(self.t + d * self.cadence())

    # ------------------------------------------------ drawing
    def draw_timeline(self):
        for m in list(STORE.months):
            if m not in self.tl_rows:
                self.tl_rows[m] = timeline_month(STORE.months[m]) + timeline_dkist(m)
        self.months_drawn = set(STORE.months)
        rows = [x for v in self.tl_rows.values() for x in v]
        df = pd.DataFrame(rows, columns=["row", "start", "end", "color", "name", "info"])
        if len(df):   # instantaneous blocks: at least a minute wide
            df["end"] = df[["end"]].assign(m=df.start + pd.Timedelta("1min")).max(axis=1)
        self.tl_src.data = {c: df[c].tolist() if c in ("row", "color", "name", "info") else df[c].values
                            for c in df.columns}

    def draw_image(self, src, fig, img, outlines, which, title, empty_note):
        if img is None:
            for k in ("img", "grid", "on", "hidden", "off", "patch", "limb", "sc_on", "sc_hidden", "mark"):
                src[k].data = {c: [] for c in src[k].data}
            fig.title.text = title
            src["note"].text = empty_note
            return
        src["note"].text = ""
        ny, nx = img["rgba"].shape[:2]
        src["img"].data = dict(image=[img["rgba"].view(np.uint32).reshape(ny, nx)], x=[img["x0"]], y=[img["y0"]],
                               dw=[img["dw"]], dh=[img["dh"]], **image_meta(img))
        gx, gy = grid(img)
        src["grid"].data = dict(xs=gx, ys=gy)
        shown = [o for o in outlines if o["show"] in self.s["show"]]
        parts = [split(o["x"], o["y"], o["kind"]) for o in shown]
        for key in KINDS:
            src[key].data = dict(xs=[p[key][0] for p in parts], ys=[p[key][1] for p in parts],
                                 color=[o["color"] for o in shown])
        ok = [o for o in shown if len(o["x"]) > 2]
        src["patch"].data = dict(xs=[o["x"] for o in ok], ys=[o["y"] for o in ok],
                                 color=[o["color"] for o in ok], label=[o["label"] for o in ok],
                                 **{f: [o["tip"][f] for o in ok] for f in ("name", "when", "cadence", "where", "soop")})
        fig.title.text = title
        r = 1.35 * img["rsun"] if which == "aia" else 1.5 * img["rsun"]   # room for off-limb polar FOVs
        for rng in (fig.x_range, fig.y_range):
            rng.reset_start, rng.reset_end = -r, r
            if self.first[which]:
                rng.start, rng.end = -r, r
        self.first[which] = False

    def draw(self):
        st = self.step
        shown_products = self.shown_products()
        blocks = st.blocks[st.blocks["product"].isin(shown_products)] if len(st.blocks) else st.blocks
        kind = "full resolution" if st.full else "synoptic"
        self.draw_image(self.src_aia, self.fig_aia, st.aia, st.outlines["aia"], "aia",
                        f"SDO/AIA {st.wl} Å  {st.aia['date'][:19].replace('T', ' ')}  ({kind}, Earth view)", "")
        has_solo = len(blocks) > 0 and st.fsi is not None
        note = ("no Solar Orbiter observation (HRI / PHI-HRT / SPICE) in this step" if not len(blocks)
                else "no FSI image within a day")
        self.draw_image(self.src_fsi, self.fig_fsi, st.fsi if has_solo else None, st.outlines["fsi"], "fsi",
                        f"SolO/EUI FSI 174 Å  {st.fsi['date'][:19].replace('T', ' ')}  (Solar Orbiter view)"
                        if has_solo else "SolO/EUI FSI 174 Å", note)
        self.draw_observers(st, blocks, has_solo)
        half = (st.t1 - st.t0) / 2
        self.title.object = (f"<b>{st.t:%Y-%m-%d %H:%M:%S} UTC</b> &nbsp; step ±{br.fmt_td(half)} &nbsp; "
                             f"cadence {st.cadence}")
        show_dk = "DKIST" in self.s["show"]
        dk_all = st.dkist if (show_dk and len(st.dkist)) else pd.DataFrame()
        lines = br.summary_lines(st, blocks, dk_all, show_dk, width=100000)
        self.summary.object = summary_html(lines)
        f = st.files
        if len(f):
            f = f[(f.start >= st.t0) & (f.start <= st.t1) & f["product"].isin(shown_products)]
            cols = [c for c in ("product", "start", "filename", "soop", "pointing_src", "study") if c in f]
            f = f[cols].sort_values("start").reset_index(drop=True)
        self.files_tab.value = f if len(f) else pd.DataFrame({"": ["no Solar Orbiter file in this step"]})
        d = dk_all
        if len(d):
            cols = [c for c in ("instrument", "start", "end", "lines", "proposal", "embargoed", "x0", "y0", "x1", "y1",
                                "frameCount", "datasetId") if c in d]
            d = d[cols].reset_index(drop=True)
        self.dkist_tab.value = d if len(d) else pd.DataFrame({"": ["no DKIST dataset in this step"]})

    def draw_observers(self, st, blocks, has_solo):
        """Sub-points and limbs: Solar Orbiter on AIA, Earth on FSI; with "positions" also PSP
        on both and the source regions of Solar Orbiter and PSP; and the positions tab."""
        def body(n, c):
            return (n, c.lat.to_value(u.deg), c.lon.wrap_at(180 * u.deg).to_value(u.deg), c.radius.to_value(u.au),
                    pf.SC_COLORS.get(n, "white"))
        have = blocks.dropna(subset=["solo_lat"]) if len(blocks) and "solo_lat" in blocks else pd.DataFrame()
        solo = None
        if len(have):
            b0 = have.iloc[0]
            solo = ("Solar Orbiter", b0.solo_lat, b0.solo_lon, b0.solo_r, pf.SC_COLORS["Solar Orbiter"])
        elif "Solar Orbiter" in st.coords:
            solo = body("Solar Orbiter", st.coords["Solar Orbiter"])
        others = [body(n, c) for n, c in st.coords.items() if n != "Solar Orbiter"]
        connect = [(n, c, pf.SC_COLORS.get(n, "white")) for n, c in st.coords.items()]
        layers = {"aia": observer_layers(st.aia, ([solo] if solo else []) + others, connect)}
        if has_solo:
            e = st.earth
            earth = ("Earth", e.lat.to_value(u.deg), 0.0, e.radius.to_value(u.au), pf.EARTH_COLOR)
            layers["fsi"] = observer_layers(st.fsi, [earth] + others, connect)
        for which, src in (("aia", self.src_aia), ("fsi", self.src_fsi)):
            if which in layers:
                for k, v in layers[which].items():
                    src[k].data = v
        if st.coords:
            data, rmax, table = positions_data(st)
            for k, v in data.items():
                self.src_pos[k].data = v
            for rng in (self.fig_pos.x_range, self.fig_pos.y_range):
                rng.start, rng.end = -rmax, rmax
                rng.reset_start, rng.reset_end = -rmax, rmax
            self.fig_pos.title.text = f"positions {st.t:%Y-%m-%d %H:%M} (Stonyhurst, au; Earth at 0°, to the right)"
            self.pos_table.object = table
            self.pos_note.visible = False
        else:
            self.pos_note.visible = True


def as_time(v):
    """A Bokeh datetime range value (ms since 1970, or a datetime) -> Timestamp."""
    return pd.Timestamp(v, unit="ms") if isinstance(v, (int, float, np.number)) else pd.Timestamp(v)


def summary_html(lines):
    """browse.summary_lines() -> HTML: each line in its colour, wrapped with a hanging indent."""
    out = []
    for text, color in lines:
        out.append(f'<div style="color:{to_hex(color)};padding-left:10ch;text-indent:-10ch;margin:1px 0">'
                   f'{html.escape(text)}</div>')
    return ('<div style="font-family:monospace;font-size:12px;white-space:pre-wrap;background:#1a1a1a;'
            'padding:8px;border-radius:4px">' + "".join(out) + "</div>")


def make_page(t=None):
    app = App(t)
    pn.state.add_periodic_callback(app.poll, period=200)
    return app.layout


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("time", nargs="?", help="UTC, e.g. 2025-10-20T12:00 (default: the last date shown)")
    ap.add_argument("--port", type=int, default=5006)
    ap.add_argument("--no-show", action="store_true", help="do not open a browser tab")
    args = ap.parse_args()
    warnings.filterwarnings("ignore", module="sunpy")
    warnings.filterwarnings("ignore", module="astropy")
    warnings.filterwarnings("ignore", category=FutureWarning)
    pn.extension("tabulator")
    pn.config.theme = "dark"
    first = [args.time]

    def page():
        t, first[0] = first[0], None   # the time given applies to the first page opened
        return make_page(t)
    pn.serve({"browse": page}, port=args.port, show=not args.no_show, title="solar_quicklook browse",
             websocket_origin=[f"localhost:{args.port}", f"127.0.0.1:{args.port}"])


if __name__ == "__main__":
    main()
