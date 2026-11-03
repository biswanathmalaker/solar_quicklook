#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
Step through time and see what was observed: an AIA image at each step, and the fields of
view of everything Solar Orbiter (HRIEUV, HRI-Lya, PHI-HRT, SPICE) and DKIST observed in
that step, with a summary of the data.

  left panel    AIA (synoptic 1024x1024, 2-min cadence; the full-resolution frame instead
                when it has been downloaded, see "original data") with the footprints of
                the Solar Orbiter FOVs (solid where AIA sees them, dotted behind the disk)
                and the DKIST FOVs
  right panel   when Solar Orbiter observed in the step: EUI/FSI 174 (nearest file within
                a day) with the FOVs as observed (dashed off the limb); otherwise empty
  summary       (bottom; each instrument in the colour of its FOV) where Solar Orbiter is;
                the campaign (SOOP) and what it is for (its ESA page); for each instrument
                what it measures, the time range of the observing block, its cadence, the
                number of files (and how many in the step), the FOV size and where it is;
                SPICE: the study, slit, raster size, exposure and the spectral lines with
                their formation temperatures (spice_meta.py); FSI: cadence of 174 and 304;
                DKIST: instrument, time, lines, FOV, proposal

A "step" is [time - cadence/2, time + cadence/2]: an observation is shown if its block
overlaps it.

  date box + Enter   go to a date (UTC: 2025-10-20, 2025-10-20 12:30, 2025-10-20T12:30:00)
  << / >>, keys b / n or left / right arrows   step back / forward by the cadence
  cadence            12s 1min 5min 30min 1h (default) 3h 6h 12h 1d 3d 1wk (AIA synoptic
                     images are every 2 min: below that, use the original data)
  wavelength         AIA channel (171 default)
  show               which observations to draw (HRIEUV, HRI-Lya, PHI-HRT, SPICE, DKIST)
  positions          a second window: spacecraft positions (polar plot, Parker spirals), and
                     AIA and FSI with PSP / Solar Orbiter / Earth sub-points, limbs and solar
                     source regions; follows the steps while open
  original data      download the full-resolution AIA frame (JSOC, 1-3 min) and the Solar
                     Orbiter files of the step (HRIEUV, HRI-Lya, PHI-HRT: the file nearest
                     to the time in each block) and show the latter in a third window; the
                     full-resolution AIA frame is then used instead of the synoptic one,
                     also next time
  toolbar            zoom / pan as usual; the zoom (in arcsec) is kept when stepping

Wavelength, cadence, the boxes ticked and the last date are kept in browse_settings.json and
restored next time. Everything downloaded or queried is kept in cache/ (see README.md):
the first step into a new month queries the Solar Orbiter archive (a few minutes, once).

Run:  ./browse.py                      # last date and settings
      ./browse.py 2025-10-20T12:00
"""
import argparse
import json
import textwrap
import threading
import traceback
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval
from matplotlib.widgets import Button, CheckButtons, RadioButtons, TextBox
from sunpy.coordinates import get_earth

import dkist_meta as dm
import plot_fov as pf
import soar_meta as sm
import soops
import spice_meta as sp

OUT_DIR = Path(__file__).resolve().parent
SETTINGS = OUT_DIR / "browse_settings.json"
AIA_FULL = OUT_DIR / "cache" / "aia_full"
SOLO_FILES = OUT_DIR / "cache" / "solo"
CADENCES = {"12s": "12s", "1min": "1min", "5min": "5min", "30min": "30min", "1h": "1h", "3h": "3h",
            "6h": "6h", "12h": "12h", "1d": "1D", "3d": "3D", "1wk": "7D"}
WAVELENGTHS = [171, 193, 211, 304, 335, 94, 131, 1600, 1700]
EUV = [171, 193, 211, 304, 335, 94, 131]   # full resolution from aia.lev1_euv_12s (download_ref_file.py)
SHOW = ["HRIEUV", "HRI-Lya", "PHI-HRT", "SPICE", "DKIST"]
PRODUCT_OF = {"HRIEUV": ["hri174", "hriopn"], "HRI-Lya": ["hrilya"], "PHI-HRT": ["phi-hrt"], "SPICE": ["spice"]}
DEFAULTS = {"wavelength": 171, "cadence": "1h", "show": SHOW, "date": "2025-10-20T12:00:00"}
FULL_TOL = pd.Timedelta("60s")   # a full-resolution AIA frame this close to the time replaces the synoptic one
FSI_WINDOW = pd.Timedelta("1D")
# what each instrument observes (summary); order of the summary
ABOUT = {
    "phi-hrt": ("PHI-HRT", "maps", "photospheric magnetic field and flows from Fe I 617.3 nm (0.5\" px): "
                "B_LOS, vector B (|B|, inclination, azimuth), v_LOS, continuum intensity"),
    "spice": ("SPICE", "files", "EUV imaging spectrograph: line intensities, Doppler shifts and widths, "
              "chromosphere to corona"),
    "hri174": ("HRIEUV", "images", "EUV images at 174 Å (Fe IX/X, corona ~1 MK), 0.49\" px"),
    "hriopn": ("HRIEUV", "images", "EUV images at 174 Å with the filter wheel open (more signal; Fe IX/X, "
               "~1 MK), 0.49\" px"),
    "hrilya": ("HRI-Lya", "images", "Lyman-alpha 121.6 nm images (H I: upper chromosphere / transition region)"),
}
FSI_ABOUT = {"fsi174": "174 Å (corona ~1 MK)", "fsi304": "304 Å (He II: chromosphere / transition region)"}
DKIST_ABOUT = {"VBI": "Visible Broadband Imager", "VISP": "Visible Spectro-Polarimeter (slit)",
               "DL-NIRSP": "Diffraction-Limited Near-IR Spectro-Polarimeter (integral field)",
               "CRYO-NIRSP": "Cryogenic Near-IR Spectro-Polarimeter (corona, limb)",
               "VTF": "Visible Tunable Filter (imaging spectropolarimeter)"}
DKIST_NAME = {"VBI": "VBI", "VISP": "ViSP", "DL-NIRSP": "DL-NIRSP", "CRYO-NIRSP": "Cryo-NIRSP", "VTF": "VTF"}
SUM_RECT = [0.005, 0.005, 0.785, 0.29]   # summary panel
SUM_BG = "0.1"


def load_settings():
    try:
        return DEFAULTS | json.loads(SETTINGS.read_text())
    except (FileNotFoundError, ValueError):
        return dict(DEFAULTS)


def save_settings(s):
    SETTINGS.write_text(json.dumps(s, indent=1))


def parse_time(text):
    return pd.Timestamp(text.strip().replace("/", "-")).tz_localize(None).floor("s")


# ---------------------------------------------------------------- metadata, by month

class MetaStore:
    """Solar Orbiter files and observing blocks (soar_meta.py) per month, kept in memory
    (the first time a month is used it may be queried from SOAR: minutes)."""

    def __init__(self):
        self.months = {}
        self.lock = threading.Lock()

    def month(self, m, status=print):
        with self.lock:
            if m in self.months:
                return self.months[m]
            end = m + pd.offsets.MonthBegin(1) - pd.Timedelta("1us")
            status(f"Solar Orbiter metadata for {m:%Y-%m} (a new month is queried from SOAR once: minutes) ...")
            f = sm.files(m, end, sm.POINTED)
            b, f = sm.blocks(f) if len(f) else (pd.DataFrame(), f)
            fsi = sm.files(m, end, ["fsi174", "fsi304"], headers=False)
            if len(f) and (f["product"] == "spice").any():
                modes = sm.spice_modes(m)
                if len(modes):
                    f = f.merge(modes, on="data_item_id", how="left")
            # block numbers are per month: make them unique
            tag = int(f"{m:%Y%m}") * 100000
            if len(b):
                b["block"] += tag
                f["block"] += tag
            self.months[m] = dict(blocks=b, files=f, fsi=fsi)
            return self.months[m]

    def window(self, t0, t1, status=print):
        """Blocks overlapping [t0, t1], their files, and the FSI files within a day."""
        ms = sm.months(t0 - FSI_WINDOW, t1 + FSI_WINDOW)
        parts = [self.month(m, status) for m in ms]
        b = pd.concat([p["blocks"] for p in parts if len(p["blocks"])], ignore_index=True) \
            if any(len(p["blocks"]) for p in parts) else pd.DataFrame()
        f = pd.concat([p["files"] for p in parts if len(p["files"])], ignore_index=True) \
            if any(len(p["files"]) for p in parts) else pd.DataFrame()
        fsi = pd.concat([p["fsi"] for p in parts if len(p["fsi"])], ignore_index=True) \
            if any(len(p["fsi"]) for p in parts) else pd.DataFrame()
        if len(b):
            b = b[(b.start <= t1) & (b.end >= t0)]
        return b, f, fsi


# ---------------------------------------------------------------- one step

def full_res_path(t, wl):
    """A downloaded full-resolution AIA frame within FULL_TOL of t, or None."""
    folder = AIA_FULL / str(wl)
    best = None
    for p in folder.glob("*.fits"):
        try:
            tp = pd.to_datetime(p.stem[-15:], format="%Y%m%d_%H%M%S")
        except ValueError:
            continue
        if abs(tp - t) <= FULL_TOL and (best is None or abs(tp - t) < abs(best[0] - t)):
            best = (tp, p)
    return None if best is None else best[1]


def load_aia(path):
    m = sunpy.map.Map(path)
    data = np.asarray(m.data, dtype=np.float32)
    if float(m.meta.get("lvl_num", 0)) == 1.0:   # full-resolution level 1: DN -> DN/s
        data = data / np.float32(m.exposure_time.to_value(u.s))
    return sunpy.map.Map(data, m.meta)


class Step:
    """Everything needed to draw one step, prepared in a background thread."""

    def __init__(self, t, cadence, wl, store, traj_names, status):
        self.t, self.cadence, self.wl = t, cadence, wl
        half = max(pd.Timedelta(CADENCES[cadence]) / 2, pd.Timedelta("1s"))
        self.t0, self.t1 = t - half, t + half
        self.blocks, self.files, fsi = store.window(self.t0, self.t1, status)
        path = full_res_path(t, wl)
        self.full = path is not None
        if path is None:
            status(f"AIA {wl} synoptic for {t:%Y-%m-%d %H:%M} ...")
            path = pf.aia_synoptic(t, wl)
        self.aia_path = path
        self.aia = load_aia(path)
        try:
            self.dkist = dm.datasets(self.t0, self.t1)
            if len(self.dkist):
                self.dkist = self.dkist.dropna(subset=["x0"])
        except Exception as e:   # network trouble: no DKIST, but the rest is shown
            status(f"DKIST: {e}")
            self.dkist = pd.DataFrame()
        self.fsi = None
        if len(self.blocks) and len(fsi):
            dt = (fsi.start - t).abs()
            k = dt.idxmin()
            if dt[k] <= FSI_WINDOW:
                status(f"FSI {fsi.start[k]:%Y-%m-%d %H:%M} ...")
                self.fsi = sunpy.map.Map(sm.download(fsi.filename[k], pf.FSI_CACHE))
        self.coords = pf.positions(t, traj_names, {}) if traj_names else {}
        self.fsi_files = fsi
        # what the campaigns are for (ESA SOOP pages) and the SPICE studies: once each, cached
        self.soop_pages, self.soop_runs, self.spice = {}, {}, {}
        for soop in dict.fromkeys(self.blocks.soop) if len(self.blocks) else []:
            if str(soop).lower() == "none":
                continue
            status(f"SOOP {soop} ...")
            self.soop_pages[soop] = soops.page(soop)
            r = soops.runs(soop, self.t0, self.t1, margin="0D")
            self.soop_runs[soop] = r.iloc[0] if len(r) else None
        if len(self.blocks) and "study" in self.files:
            sp_files = self.files[self.files.block.isin(self.blocks.block[self.blocks["product"] == "spice"])]
            for study, g in sp_files.dropna(subset=["study"]).groupby("study"):
                in_step = g[(g.start >= self.t0) & (g.start <= self.t1)]
                if len(in_step) or (g.start - t).abs().min() <= pd.Timedelta(CADENCES[cadence]):
                    status(f"SPICE study {study} ...")
                    self.spice[study] = sp.study(study, g.filename.iloc[0])


# ---------------------------------------------------------------- the browser

class Browser:
    def __init__(self, t=None):
        self.s = load_settings()
        self.t = parse_time(t or self.s["date"])
        self.store = MetaStore()
        self.pool = ThreadPoolExecutor(max_workers=1)    # steps (and their prefetch)
        self.slow = ThreadPoolExecutor(max_workers=1)    # original data downloads
        self.pending = None      # (key, future) of the step being prepared
        self.prefetched = {}     # key -> future
        self.step = None
        self.status_text = ""
        self.view = {}           # panel -> limits in arcsec, kept when stepping
        self.quiet = False       # True while go() sets the date box
        self.prev_aia = self.prev_fsi = None
        self.pos_fig = None
        self.orig_fig = None
        self.orig_job = None

        self.fig = plt.figure(figsize=(16, 9.2))
        self.fig.canvas.manager.set_window_title("solar_quicklook: browse")
        self.rect_aia = [0.04, 0.365, 0.33, 0.56]
        self.rect_fsi = [0.43, 0.365, 0.32, 0.56]
        self.ax_aia = self.ax_fsi = None
        self.title = self.fig.text(0.40, 0.99, "", ha="center", va="top", fontsize=12)
        self.ax_sum = self.fig.add_axes(SUM_RECT)
        self.status = self.fig.text(0.40, 0.305, "", fontsize=9, color="tab:red", ha="center")

        # controls
        self.date_box = TextBox(self.fig.add_axes([0.80, 0.905, 0.17, 0.04]), "date ", initial=f"{self.t:%Y-%m-%d %H:%M:%S}")
        self.date_box.on_submit(self.date_typed)
        self.back = Button(self.fig.add_axes([0.80, 0.85, 0.08, 0.04]), "<<")
        self.fwd = Button(self.fig.add_axes([0.89, 0.85, 0.08, 0.04]), ">>")
        self.back.on_clicked(lambda _: self.go(self.t - pd.Timedelta(CADENCES[self.s["cadence"]])))
        self.fwd.on_clicked(lambda _: self.go(self.t + pd.Timedelta(CADENCES[self.s["cadence"]])))
        ax = self.fig.add_axes([0.80, 0.44, 0.08, 0.39], title="cadence")
        self.cad = RadioButtons(ax, list(CADENCES), active=list(CADENCES).index(self.s["cadence"]))
        self.cad.on_clicked(self.cadence_set)
        ax = self.fig.add_axes([0.89, 0.50, 0.08, 0.33], title="AIA (Å)")
        wls = [str(w) for w in WAVELENGTHS]
        self.wl = RadioButtons(ax, wls, active=wls.index(str(self.s["wavelength"])))
        self.wl.on_clicked(self.wavelength_set)
        for lab in self.cad.labels + self.wl.labels:
            lab.set_fontsize(8)
        ax = self.fig.add_axes([0.80, 0.25, 0.17, 0.16], title="show")
        self.show = CheckButtons(ax, SHOW, [n in self.s["show"] for n in SHOW])
        self.show.on_clicked(self.show_set)
        self.pos_button = Button(self.fig.add_axes([0.80, 0.18, 0.17, 0.045]), "positions")
        self.pos_button.on_clicked(self.open_positions)
        self.orig_button = Button(self.fig.add_axes([0.80, 0.12, 0.17, 0.045]), "original data")
        self.orig_button.on_clicked(self.original_data)

        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.fig.canvas.mpl_connect("close_event", self.on_close)
        self.timer = self.fig.canvas.new_timer(interval=250)
        self.timer.add_callback(self.poll)
        self.timer.start()
        self.go(self.t)

    # ------------------------------------------------ settings and controls
    def set_status(self, text):
        self.status_text = text   # from worker threads: drawn by poll()

    def date_typed(self, text):
        if self.quiet:   # set by go(), not typed
            return
        try:
            self.go(parse_time(text))
        except ValueError:
            self.set_status(f"can't read the date '{text}'")

    def cadence_set(self, label):
        self.s["cadence"] = label
        save_settings(self.s)
        self.prefetched.clear()
        self.go(self.t)

    def wavelength_set(self, label):
        self.s["wavelength"] = int(label)
        save_settings(self.s)
        self.prefetched.clear()
        self.go(self.t)

    def show_set(self, _):
        self.s["show"] = [n for n, on in zip(SHOW, self.show.get_status()) if on]
        save_settings(self.s)
        if self.step is not None:
            self.draw()

    def on_key(self, event):
        if event.key in ("b", "left"):
            self.go(self.t - pd.Timedelta(CADENCES[self.s["cadence"]]))
        elif event.key in ("n", "right"):
            self.go(self.t + pd.Timedelta(CADENCES[self.s["cadence"]]))

    def on_close(self, _):
        self.s["date"] = f"{self.t:%Y-%m-%dT%H:%M:%S}"
        save_settings(self.s)
        self.timer.stop()
        self.pool.shutdown(wait=False, cancel_futures=True)
        self.slow.shutdown(wait=False, cancel_futures=True)
        for f in (self.pos_fig, self.orig_fig):
            if f is not None:
                plt.close(f)

    # ------------------------------------------------ stepping
    def traj_names(self):
        return ["Solar Orbiter", "PSP"] if self.pos_fig is not None else []

    def key(self, t):
        return (t, self.s["cadence"], self.s["wavelength"], bool(self.traj_names()))

    def submit(self, t):
        k = self.key(t)
        if k not in self.prefetched:
            def status(text):   # only the step being waited for reports (not a prefetch)
                if self.pending is not None and self.pending[0] == k:
                    self.set_status(text)
            self.prefetched[k] = self.pool.submit(Step, t, self.s["cadence"], self.s["wavelength"],
                                                  self.store, self.traj_names(), status)
        return k, self.prefetched[k]

    def go(self, t):
        self.t = t
        self.quiet = True
        self.date_box.set_val(f"{t:%Y-%m-%d %H:%M:%S}")
        self.quiet = False
        self.pending = self.submit(t)
        self.set_status(f"loading {t:%Y-%m-%d %H:%M:%S} ...")
        self.poll()

    def poll(self):
        if self.status.get_text() != self.status_text:
            self.status.set_text(self.status_text)
            self.fig.canvas.draw_idle()
        if self.pending is not None and self.pending[1].done():
            k, fut = self.pending
            self.pending = None
            try:
                self.step = fut.result()
            except Exception as e:
                traceback.print_exc()
                self.set_status(f"failed: {e}")
                self.prefetched.pop(k, None)
                return
            self.set_status("")
            self.draw()
            # keep the cache of prepared steps small, then prefetch the next one
            for old in [x for x in self.prefetched if x[0] not in (self.t, self.t + pd.Timedelta(CADENCES[self.s["cadence"]]))]:
                self.prefetched.pop(old)
            self.submit(self.t + pd.Timedelta(CADENCES[self.s["cadence"]]))
        if self.orig_job is not None and self.orig_job.done():
            job, self.orig_job = self.orig_job, None
            try:
                t, solo_paths = job.result()
            except Exception as e:
                traceback.print_exc()
                self.set_status(f"original data failed: {e}")
                return
            self.set_status("")
            self.prefetched.clear()
            if t == self.t:
                self.go(self.t)            # now with the full-resolution AIA frame
            self.show_originals(t, solo_paths)

    # ------------------------------------------------ drawing
    def remember_view(self, name, ax, m):
        """Keep the zoom of a panel (in arcsec) unless it shows the whole disk."""
        if ax is None or m is None:
            return
        x, y = ax.get_xlim(), ax.get_ylim()
        c = m.wcs.pixel_to_world([x[0], x[1]], [y[0], y[1]])
        r = m.rsun_obs.to_value(u.arcsec)
        tx, ty = c.Tx.to_value(u.arcsec), c.Ty.to_value(u.arcsec)
        whole = np.ptp(tx) > 2.2 * r
        self.view[name] = None if whole else (tx, ty)

    def restore_view(self, name, ax, m):
        v = self.view.get(name)
        if not v:
            return
        px = m.wcs.world_to_pixel(SkyCoord(v[0] * u.arcsec, v[1] * u.arcsec, frame=m.coordinate_frame))
        ax.set_xlim(*px[0])
        ax.set_ylim(*px[1])

    def new_axes(self, rect, m, old):
        if old is not None:
            old.remove()
        return self.fig.add_axes(rect, projection=m) if m is not None else self.fig.add_axes(rect)

    def shown_blocks(self):
        st = self.step
        prods = [p for n in self.s["show"] for p in PRODUCT_OF.get(n, [])]
        return st.blocks[st.blocks["product"].isin(prods)] if len(st.blocks) else st.blocks

    def draw(self):
        st = self.step
        if self.ax_aia is not None:
            self.remember_view("aia", self.ax_aia, self.prev_aia)
            self.remember_view("fsi", self.ax_fsi, self.prev_fsi)
        self.ax_aia = self.new_axes(self.rect_aia, st.aia, self.ax_aia)
        pf.show_map(self.ax_aia, st.aia)
        kind = "full resolution" if st.full else "synoptic"
        self.ax_aia.set_title(f"SDO/AIA {st.wl} Å  {st.aia.date.strftime('%Y-%m-%d %H:%M:%S')}  ({kind})", fontsize=10)
        blocks = self.shown_blocks()
        items = []
        for _, b in blocks.iterrows():
            if np.isfinite(b.get("c0x", np.nan)):
                pf.draw_fov_earth(self.ax_aia, st.aia, b, pf.COLORS[b["product"]])
        for p in dict.fromkeys(blocks["product"]) if len(blocks) else []:
            if (pf.LABELS[p], pf.COLORS[p], "-") not in items:
                items.append((pf.LABELS[p], pf.COLORS[p], "-"))
        dk = st.dkist if ("DKIST" in self.s["show"] and len(st.dkist)) else pd.DataFrame()
        if len(dk):
            key = dk[["x0", "y0", "x1", "y1"]].round(-1).astype(str).agg(",".join, axis=1)
            dk = dk.assign(key=dk.instrument + key).drop_duplicates("key")
            items += [(f"DKIST {i}", pf.DKIST_COLORS.get(str(i).upper(), pf.DKIST_OTHER), "-") for i in dict.fromkeys(dk.instrument)]

        has_solo = len(blocks) > 0 and st.fsi is not None
        self.ax_fsi = self.new_axes(self.rect_fsi, st.fsi if has_solo else None, self.ax_fsi)
        if has_solo:
            pf.show_map(self.ax_fsi, st.fsi, extent=1.5)   # room for off-limb polar FOVs
            self.ax_fsi.set_title(f"SolO/EUI FSI 174 Å  {st.fsi.date.strftime('%Y-%m-%d %H:%M:%S')}", fontsize=10)
            for _, b in blocks.iterrows():
                if np.isfinite(b.get("c0x", np.nan)):
                    corners = b[sm.CORNER_COLS].to_numpy(float).reshape(4, 2)
                    tx, ty = pf.box_outline(corners)
                    pf.draw_hpc_box(self.ax_fsi, st.fsi, tx, ty, b.solo_lat, b.solo_lon, b.solo_r, pf.COLORS[b["product"]])
            pf.mark_observer(self.ax_fsi, st.fsi, get_earth(st.fsi.date), pf.EARTH_COLOR, "Earth")
            self.restore_view("fsi", self.ax_fsi, st.fsi)
        else:
            self.ax_fsi.set_axis_off()
            self.ax_fsi.text(0.5, 0.5, "no Solar Orbiter observation (HRI / PHI-HRT / SPICE) in this step"
                             if not len(blocks) else "no FSI image within a day",
                             ha="center", va="center", transform=self.ax_fsi.transAxes, color="0.4")
        if len(dk):
            pf.draw_dkist(self.ax_aia, st.aia, self.ax_fsi if has_solo else None, st.fsi, dk)
        if blocks is not None and len(blocks) and np.isfinite(blocks.solo_lat).any():
            b0 = blocks.dropna(subset=["solo_lat"]).iloc[0]
            pf.mark_observer(self.ax_aia, st.aia, pf.solo_observer(b0),
                             pf.SC_COLORS["Solar Orbiter"], "Solar Orbiter")
        pf.legend(self.ax_aia, items)
        if has_solo:
            pf.legend(self.ax_fsi, items)
        self.restore_view("aia", self.ax_aia, st.aia)
        self.prev_aia, self.prev_fsi = st.aia, (st.fsi if has_solo else None)

        half = (st.t1 - st.t0) / 2
        self.title.set_text(f"{st.t:%Y-%m-%d %H:%M:%S} UTC   step ±{fmt_td(half)}   cadence {st.cadence}")
        show_dk = "DKIST" in self.s["show"]
        dk_all = st.dkist if (show_dk and len(st.dkist)) else pd.DataFrame()
        self.draw_summary(lambda width: summary_lines(st, blocks, dk_all, show_dk, width))
        self.fig.canvas.draw_idle()
        if self.pos_fig is not None:
            self.draw_positions()

    def draw_summary(self, make):
        """make(width) -> [(text, colour)] wrapped to width characters; drawn top down in the
        summary panel with the largest font (<= 8 pt) at which they fit."""
        ax = self.ax_sum
        ax.cla()
        ax.set_facecolor(SUM_BG)
        ax.set_xticks([])
        ax.set_yticks([])
        w_pt = SUM_RECT[2] * self.fig.get_figwidth() * 72
        h_pt = SUM_RECT[3] * self.fig.get_figheight() * 72 - 4
        for fs in np.arange(8.0, 5.4, -0.25):
            lines = make(int(w_pt / (0.602 * fs)) - 3)
            if 1.25 * fs * len(lines) <= h_pt:
                break
        dy = 1.25 * fs / (h_pt + 4)
        for k, (text, color) in enumerate(lines):
            ax.text(0.003, 1 - 2 / (h_pt + 4) - (k + 0.5) * dy, text, color=color, family="monospace",
                    fontsize=fs, va="center", ha="left", transform=ax.transAxes, clip_on=True)

    # ------------------------------------------------ positions window
    def open_positions(self, _):
        if self.pos_fig is not None and plt.fignum_exists(self.pos_fig.number):
            self.pos_fig.canvas.manager.show()
            return
        self.pos_fig = plt.figure(figsize=(17, 6.8))
        self.pos_fig.canvas.manager.set_window_title("solar_quicklook: positions")
        self.pos_fig.canvas.mpl_connect("close_event", self.positions_closed)
        self.prefetched.clear()
        self.go(self.t)     # the step again, now with the spacecraft positions
        plt.show(block=False)

    def positions_closed(self, _):
        self.pos_fig = None

    def draw_positions(self):
        st, fig = self.step, self.pos_fig
        if not st.coords:
            return
        fig.clf()
        gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.8])
        v_range = pf.V_RANGE
        ax1 = fig.add_subplot(gs[0], projection=st.aia)
        pf.show_map(ax1, st.aia)
        ax1.set_title(f"AIA {st.wl} Å (Earth view)", fontsize=10)
        ax2 = None
        if st.fsi is not None:
            ax2 = fig.add_subplot(gs[1], projection=st.fsi)
            pf.show_map(ax2, st.fsi)
            ax2.set_title(f"FSI 174 Å {st.fsi.date.strftime('%Y-%m-%d %H:%M')} (Solar Orbiter view)", fontsize=10)
            pf.mark_observer(ax2, st.fsi, get_earth(st.fsi.date), pf.EARTH_COLOR, "Earth")
        items = []
        for n, c in st.coords.items():
            if c is None:
                continue
            color = pf.SC_COLORS.get(n, "white")
            pf.mark_observer(ax1, st.aia, c, color, n)
            if ax2 is not None and n != "Solar Orbiter":
                pf.mark_observer(ax2, st.fsi, c, color, n)
            for ax, m in ((ax1, st.aia), (ax2, st.fsi)):
                if ax is not None:
                    pf.draw_connectivity(ax, m, c, n, color, v_range, pf.MARGIN)
            items.append((f"{n} source region", color, "-"))
        pf.legend(ax1, items)
        pf.draw_positions(fig.add_subplot(gs[2], projection="polar"), st.t, st.coords, list(st.coords), v_range)
        fig.suptitle(f"{st.t:%Y-%m-%d %H:%M} UTC   source regions: Parker spiral to 2.5 Rsun, "
                     f"{v_range[0]}-{v_range[1]} km/s ±{pf.MARGIN}°, no PFSS", fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.canvas.draw_idle()

    # ------------------------------------------------ original data
    def original_data(self, _):
        if self.orig_job is not None:
            self.set_status("original data: already downloading ...")
            return
        st = self.step
        if st is None:
            return
        blocks = self.shown_blocks()
        files = st.files
        wanted = []
        for _, b in blocks.iterrows():
            if b["product"] == "spice":
                continue   # SPICE rasters are large: not fetched here
            g = files[files.block == b.block]
            if len(g):
                wanted.append(g.iloc[(g.start - st.t).abs().argmin()])
        self.set_status(f"original data: AIA {st.wl} full resolution (JSOC, 1-3 min) + {len(wanted)} Solar Orbiter file(s) ...")
        self.orig_job = self.slow.submit(fetch_originals, st.t, st.wl, st.full, wanted, self.set_status)

    def show_originals(self, t, paths):
        if not paths:
            return
        if self.orig_fig is None or not plt.fignum_exists(self.orig_fig.number):
            self.orig_fig = plt.figure(figsize=(6.5 * len(paths), 6.8))
            self.orig_fig.canvas.manager.set_window_title("solar_quicklook: Solar Orbiter data")
        fig = self.orig_fig
        fig.clf()
        fig.set_size_inches(6.5 * len(paths), 6.8)
        for k, (product, path) in enumerate(paths):
            m = sunpy.map.Map(path)
            ax = fig.add_subplot(1, len(paths), k + 1, projection=m)
            data = np.where(np.isfinite(m.data), m.data, 0)
            if product == "phi-hrt":
                m.plot(axes=ax, cmap="gray", vmin=-100, vmax=100, title=False)
            else:
                norm = ImageNormalize(data, interval=PercentileInterval(99.7), stretch=AsinhStretch(0.05))
                m.plot(axes=ax, norm=norm, title=False)
            m.draw_grid(axes=ax, grid_spacing=5 * u.deg, color="white", alpha=0.3, lw=0.5)
            ax.set_title(f"{pf.LABELS[product]}  {m.date.strftime('%Y-%m-%d %H:%M:%S')}\n{Path(path).name}", fontsize=9)
        fig.suptitle(f"Solar Orbiter files nearest to {t:%Y-%m-%d %H:%M:%S}"
                     f"{'  (PHI-HRT blos, ±100 G)' if any(p == 'phi-hrt' for p, _ in paths) else ''}", fontsize=10)
        fig.tight_layout()
        fig.canvas.draw_idle()
        plt.show(block=False)


# ---------------------------------------------------------------- the summary

WIDTH = 190   # characters per summary line
DIM = "0.8"   # colour of the continuation lines


def wrap(text, first, indent=10, width=WIDTH, max_lines=None):
    """[(line, colour)]: text wrapped, the first line in colour first, the rest dimmer."""
    out = textwrap.wrap(text, width=width, subsequent_indent=" " * indent)
    if max_lines and len(out) > max_lines:
        out = out[:max_lines]
        out[-1] = out[-1][:width - 3] + "..."
    return [(line, first if k == 0 else DIM) for k, line in enumerate(out)]


def west_east(lon):
    return f"{abs(lon):.1f}° {'west' if lon >= 0 else 'east'} of the Sun-Earth line"


def fov_text(b):
    """Size of the FOV (arcsec as seen from Solar Orbiter, and Mm on the Sun) and where it is."""
    if not np.isfinite(b.get("c0x", np.nan)):
        return "pointing unknown"
    c = np.array([[b[f"c{k}x"], b[f"c{k}y"]] for k in range(4)], float)
    w, h = np.hypot(*(c[1] - c[0])), np.hypot(*(c[3] - c[0]))
    mm = np.radians(1 / 3600) * b.solo_r * sm.AU_KM / 1e3
    out = f"FOV {w:.0f}\"x{h:.0f}\" ({w * mm:.0f}x{h * mm:.0f} Mm)"
    if np.isfinite(b.get("lat", np.nan)):
        out += f", centre lat {b.lat:+.0f}° lon {b.lon:+.0f}°"
        if np.isfinite(b.get("lat_min", np.nan)):
            out += f", covers lat {b.lat_min:+.0f}..{b.lat_max:+.0f}°"
        if np.isfinite(b.get("earth_vis", np.nan)):
            out += f"; {100 * b.earth_vis:.0f}% also seen from Earth"
    else:
        out += ", centre off the disk"
    return out


def spice_text(st, files):
    """Studies of the SPICE files (those in the step, else all): description, raster, lines."""
    out = []
    if "study" not in files:
        return out
    for study, g in files.dropna(subset=["study"]).groupby("study"):
        info = st.spice.get(study)
        if info is None:
            if study in st.spice:   # tried, not readable
                out.append(f"study {study} ({len(g)} files; windows unknown: not in a SPICE data release)")
            continue
        kind = str(info.get("studytyp") or "")
        nx, slit, cd1 = info.get("nx") or 1, info.get("slit_wid"), info.get("cdelt1") or 0
        ny, cd2 = info.get("ny") or 0, info.get("cdelt2") or 0
        geom = (f"raster of {nx} slit positions x {cd1:g}\" = {nx * cd1:.0f}\" by {ny * cd2:.0f}\" along the slit"
                if nx > 1 else f"sit-and-stare, {ny * cd2:.0f}\" along the slit")
        out.append(f"study \"{info.get('studydes') or study}\" [{study}]: {kind.lower()}, {slit}\" slit, {geom}, "
                   f"{info.get('xposure')} s exposures, {fmt_s(float(info.get('telapse') or np.nan))} per file"
                   f"{'; ' + info['obs_desc'] if info.get('obs_desc') else ''}")
        out.append(f"lines (formation T): {sp.lines_text(info)}")
    return out


def fsi_text(st):
    """FSI 174 and 304 around the time: cadence, number in the step, or the nearest image."""
    f = st.fsi_files
    if f is None or not len(f):
        return None
    parts = []
    for prod, about in FSI_ABOUT.items():
        g = f[f["product"] == prod].sort_values("start")
        if not len(g):
            continue
        near = g[(g.start - st.t).abs() <= pd.Timedelta("12h")]
        cad = near.start.diff().dt.total_seconds().median() if len(near) > 1 else np.nan
        n_in = int(((g.start >= st.t0) & (g.start <= st.t1)).sum())
        if n_in:
            parts.append(f"{about}: one every {fmt_s(cad)}, {n_in} in this step")
        else:
            k = (g.start - st.t).abs().idxmin()
            parts.append(f"{about}: none in this step, nearest {g.start[k]:%m-%d %H:%M}")
    return "full-disk EUV images, " + "; ".join(parts) if parts else None


def summary_lines(st, blocks, dk, show_dkist, width=WIDTH):
    """The summary of a step: [(text, colour)], lines of up to width characters."""
    def wr(text, first, **kw):
        return wrap(text, first, width=width, **kw)
    lines = []
    solo = pf.SC_COLORS["Solar Orbiter"]
    have = blocks.dropna(subset=["solo_lat"]) if len(blocks) else blocks
    step = f"this step {st.t0:%m-%d %H:%M:%S} - {st.t1:%m-%d %H:%M:%S} UTC"
    if len(have):
        b0 = have.iloc[0]
        lines.append((f"SOLAR ORBITER  {b0.solo_r:.2f} au from the Sun, latitude {b0.solo_lat:+.1f}°, "
                      f"{west_east(b0.solo_lon)}        ({step})", solo))
    else:
        lines.append((f"SOLAR ORBITER  no HRIEUV / HRI-Lya / PHI-HRT / SPICE observation in {step}", solo))
    # campaigns
    for soop in dict.fromkeys(blocks.soop) if len(blocks) else []:
        if str(soop).lower() == "none":
            continue
        run, page = st.soop_runs.get(soop), st.soop_pages.get(soop)
        head = f"Campaign (SOOP) {str(soop).split('_')[-1]}  [{soop}]"
        if run is not None:
            head += f"  run {run.start:%Y-%m-%d} - {run.end:%Y-%m-%d}"
            head += f", {run.rsw}" if str(run.rsw).strip() else ""
            head += f", coordinators: {run.coordinators}" if str(run.coordinators).strip() else ""
            coord = [f"{k} {run[c]}" for k, c in (("IRIS", "iris"), ("Hinode", "hinode"), ("DKIST", "dkist"))
                     if str(run[c]).strip() not in ("", "nan")]
            head += f"; coordinated with {'; '.join(coord)}" if coord else ""
        lines.append((head, solo))
        if page and page.get("description"):
            lines += wr("aim: " + page["description"], DIM, indent=5, max_lines=2)
    # instruments
    n_in = st.files[(st.files.start >= st.t0) & (st.files.start <= st.t1)] if len(st.files) else st.files
    for prod, (name, unit, about) in ABOUT.items():
        bs = blocks[blocks["product"] == prod].sort_values("start") if len(blocks) else blocks
        if not len(bs):
            continue
        color = pf.COLORS[prod]
        lines += wr(f"{name:9} {about}", color)
        for k, (_, b) in enumerate(bs.iterrows()):
            if k == 2:
                lines.append((f"{'':10}+ {len(bs) - 2} more {name} block(s) in this step", DIM))
                break
            k_in = int((n_in.block == b.block).sum()) if len(n_in) else 0
            u = "rasters" if prod == "spice" and set(str(b.descriptors).split(",")) <= {"spice-n-ras", "spice-w-ras"} else unit
            cad = f"one every {fmt_s(b.cadence_s)}" if int(b.n) > 1 else "a single file"
            lines += wr(f"{'':10}{b.start:%m-%d %H:%M} - {b.end:%m-%d %H:%M} ({fmt_td(b.end - b.start)}): "
                          f"{cad}, {int(b.n)} {u}, {k_in} in this step; {fov_text(b)}"
                          f"{'; no SOOP name in the archive' if str(b.soop).lower() == 'none' else ''}", DIM)
            if prod == "spice":
                g = st.files[st.files.block == b.block]
                g_in = g[(g.start >= st.t0) & (g.start <= st.t1)]
                for text in spice_text(st, g_in if len(g_in) else g):
                    lines += wr(f"{'':10}{text}", DIM)
    fsi = fsi_text(st)
    if fsi:
        lines += wr(f"{'FSI':9} {fsi}", "0.95")
    # DKIST: one line per instrument and proposal
    if show_dkist:
        if len(dk):
            for k, ((inst, prop), g) in enumerate(dk.groupby([dk.instrument.str.upper(), "proposal"], sort=False)):
                if k == 4:
                    lines.append((f"{'':10}+ more DKIST instruments / proposals", DIM))
                    break
                color = pf.DKIST_COLORS.get(inst, pf.DKIST_OTHER)
                x0, y0, x1, y1 = g.x0.min(), g.y0.min(), g.x1.max(), g.y1.max()
                t0, t1 = g.start.min(), g.end.max()
                when = f"{t0:%m-%d %H:%M} - {t1:%H:%M}" if t0.date() == t1.date() else f"{t0:%m-%d %H:%M} - {t1:%m-%d %H:%M}"
                lns = ", ".join(dict.fromkeys(str(x) for x in g.lines.dropna()))
                extra = []
                if "frameCount" in g and g.frameCount.notna().any():
                    extra.append(f"{int(g.frameCount.sum())} frames")
                if "averageDatasetSpatialSampling" in g and g.averageDatasetSpatialSampling.notna().any():
                    extra.append(", ".join(f"{v:.3g}" for v in sorted(set(g.averageDatasetSpatialSampling.dropna().round(4))))
                                 + '" px')
                emb = ""
                if g.embargoed.any():
                    until = str(g.embargoEndDate.dropna().max())[:10] if "embargoEndDate" in g else ""
                    emb = f", embargoed{' until ' + until if until else ''}"
                lines += wr(f"{'DKIST':9} {DKIST_NAME.get(inst, inst)}, {DKIST_ABOUT.get(inst, '')}: {len(g)} dataset(s) "
                            f"{when}; lines: {lns}; FOV up to {x1 - x0:.0f}\"x{y1 - y0:.0f}\" at ({(x0 + x1) / 2:+.0f}\", "
                            f"{(y0 + y1) / 2:+.0f}\") as seen from Earth; {'; '.join(extra)}; proposal {prop}{emb}", color)
        else:
            lines.append(("DKIST     no dataset in this step", DIM))
    return lines


def fetch_originals(t, wl, have_full, rows, status):
    """Full-resolution AIA frame (unless there) and the Solar Orbiter files; returns
    (t, [(product, path)])."""
    import download_ref_file as drf
    if not have_full and wl not in EUV:
        status(f"original data: AIA {wl} full resolution not handled (UV series); Solar Orbiter files only")
    elif not have_full:
        status(f"original data: AIA {wl} full resolution from JSOC (1-3 min) ...")
        try:
            drf.download("aia", wl, Time(t.to_pydatetime()), folder=AIA_FULL / str(wl))
        except Exception as e:
            status(f"AIA full resolution failed: {e}")
    paths = []
    for r in rows:
        status(f"original data: {r.filename} ...")
        paths.append((r["product"], sm.download(r.filename, SOLO_FILES / r["product"])))
    return t, paths


def fmt_td(td):
    s = td.total_seconds()
    return f"{s:.0f}s" if s < 120 else f"{s / 60:.0f}min" if s < 7200 else f"{s / 3600:.1f}h" if s < 172800 else f"{s / 86400:.1f}d"


def fmt_s(s):
    return "-" if not np.isfinite(s) else fmt_td(pd.Timedelta(seconds=s))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("time", nargs="?", help="UTC, e.g. 2025-10-20T12:00 (default: the last date shown)")
    args = ap.parse_args()
    warnings.filterwarnings("ignore", module="sunpy")
    warnings.filterwarnings("ignore", module="astropy")
    warnings.filterwarnings("ignore", category=FutureWarning)
    Browser(args.time)
    plt.show()


if __name__ == "__main__":
    main()
