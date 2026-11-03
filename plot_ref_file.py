#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
Show the AIA frames of {name}/{wavelength}/ (from download_ref_file.py: one file, or
several to step through, in time order) and pick the plume's reference coordinate
by clicking on it.

  left click     mark a point; its helioprojective coordinate (arcsec, at the
                 file's DATE-OBS) is printed and written on the figure
  cursor         the toolbar shows helioprojective x, y (arcsec) and heliographic
                 Stonyhurst lon, lat (deg) under the mouse
  u              remove the last point
  sliders        vmin / vmax of the display (DN/s); the image is redrawn when
                 the mouse button is released
  stretch        linear (default), sqrt, log, asinh or squared colour scale
  inset          then click two opposite corners of a region on the image (or type
                 "x, y" in arcsec in the bottom left / top right boxes + Enter): a
                 zoomed inset of it appears, joined to the region by lines; drag
                 it to move it, drag a corner to resize it; remove inset removes it
  spacecraft     tick PSP / Solar Orbiter / STEREO-A / BepiColombo (right panel) to
                 draw its sub-spacecraft point (the centre of the disk it sees;
                 open marker at the limb if on the far side) and its limb (solid
                 where AIA sees it, dotted where hidden); positions from JPL
                 Horizons at DATE-OBS, as in plot_positions.py. Earth / L1 (Wind,
                 ACE) has only the connectivity region
  connectivity   (tick box + speed range + margin, right panel) for each ticked
                 spacecraft, the region of the Sun it can be magnetically
                 connected to: its latitude, and the longitudes a Parker spiral
                 for the speed range (default 300-800 km/s) reaches at 2.5 Rsun
                 (west of the sub-spacecraft point), widened by the margin
                 (default 10 deg), drawn radially on the disk. No PFSS: the
                 photospheric footpoint can lie elsewhere in a coronal hole.
                 The style button cycles: outline (default) / shaded / dim
                 outside (the image outside all regions darkened)
  prev / next    (with several frames) step through them; also keys b / n. The
                 image, grid, points, inset and spacecraft views follow;
                 contrast and zoom are kept
  save png/pdf   save the image (with the marked points, without the widgets)
                 next to the frame's .fits file, with the same name
  toolbar        zoom / pan as usual (clicks are ignored while zoom or pan is on)

When the window is closed, the last point is saved to ref_coords.csv in this
folder (one row per plume: plume_name, date_obs, xcen, ycen), replacing the
plume's earlier row if there is one.

Run:  ./plot_ref_file.py --name plume_002              # plume_002/171/
      ./plot_ref_file.py --name plume_002 --wavelength 193
"""
import argparse
import csv
import gc
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from matplotlib.colors import to_rgba
from matplotlib.path import Path as MplPath
from matplotlib.patches import Polygon
from sunpy.coordinates import HeliographicStonyhurst, get_body_heliographic_stonyhurst
from sunpy.coordinates.sun import L0 as carrington_L0
from astropy.visualization import (AsinhStretch, ImageNormalize, LinearStretch, LogStretch,
                                   SqrtStretch, SquaredStretch)
from matplotlib import patheffects
from matplotlib.backend_tools import Cursors
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.widgets import Button, CheckButtons, RadioButtons, Slider, TextBox
from sunpy.visualization.drawing import limb

from plot_positions import R_SS, V_RANGE, Trajectories, parse_speeds, spiral_shift

OUT_DIR = Path(__file__).resolve().parent
COORDS_CSV = OUT_DIR / "ref_coords.csv"
FIELDS = ["plume_name", "date_obs", "xcen", "ycen"]
WAVELENGTHS = [171, 193, 211, 304, 335, 94, 131]
STRETCHES = {   # first one is the default
    "linear": LinearStretch,
    "sqrt": SqrtStretch,
    "log": LogStretch,
    "asinh": lambda: AsinhStretch(0.02),
    "squared": SquaredStretch,
}
INSET_COLOR = "white"
BUTTON_COLOR, ACTIVE_COLOR = "0.85", "gold"
# colours that stand out on the AIA colour maps and differ from the picked points (red / cyan)
VIEWS = {"PSP": "magenta", "Solar Orbiter": "lime", "STEREO-A": "deepskyblue", "BepiColombo": "blue",
         "Earth / L1": "red"}
L1 = "Earth / L1"   # Wind, ACE: from the ephemeris, at 0.99 of Earth's distance; connectivity region only
MARGIN = 10         # deg, default margin around the connectivity region
REGION_STYLES = ["outline", "shaded", "dim outside"]   # the first is the default
DIM = 0.6           # "dim outside": opacity of the black over the image outside the regions
VIEWS_SHOWN = ["PSP", "Solar Orbiter"]   # ticked at the start
OUTLINE = [patheffects.withStroke(linewidth=2, foreground="black")]
AX_POS = [0.08, 0.21, 0.7, 0.71]
HANDLE_PX = 15        # grab an inset corner within this distance (pixels) to resize it
MIN_INSET_PX = 40


def convex_hull(points):
    """Andrew's monotone chain; points are (x, y, tag) tuples; returns the hull in order."""
    pts = sorted(set(points), key=lambda q: (q[0], q[1]))

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for q in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], q) <= 0:
            lower.pop()
        lower.append(q)
    for q in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], q) <= 0:
            upper.pop()
        upper.append(q)
    return lower[:-1] + upper[:-1]


def data_dir(name, wavelength):
    return OUT_DIR / name / str(wavelength)


def frame_time(path):
    """DATE-OBS from the file name ({name}_{wavelength}_YYYYMMDD_HHMMSS.fits), or the header."""
    match = re.search(r"_(\d{8}_\d{6})\.fits$", path.name)
    if match:
        return match.group(1)
    with fits.open(path) as hdul:
        return Time(hdul[-1].header["DATE-OBS"]).strftime("%Y%m%d_%H%M%S")


def load_map(path):
    m = sunpy.map.Map(path)
    # DN/s; float32 (64 MB per frame instead of 128 MB) is plenty for AIA
    data = np.asarray(m.data, dtype=np.float32)   # no copy: AIA level 1 is already float32
    return sunpy.map.Map(data / np.float32(m.exposure_time.to_value(u.s)), m.meta)


def save_coord(plume_name, date_obs, x, y):
    rows = []
    if COORDS_CSV.is_file():
        with open(COORDS_CSV, newline="") as f:
            rows = [r for r in csv.DictReader(f) if r["plume_name"] != plume_name]
    rows.append({"plume_name": plume_name, "date_obs": date_obs, "xcen": f"{x:.1f}", "ycen": f"{y:.1f}"})
    rows.sort(key=lambda r: r["plume_name"])
    tmp = COORDS_CSV.with_suffix(".tmp")
    with open(tmp, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(COORDS_CSV)
    print(f"Saved {plume_name}: xcen={x:.1f}, ycen={y:.1f} -> {COORDS_CSV}")


class QuickTextBox(TextBox):
    """TextBox that repaints only its own box (blitting): the stock one redraws the
    whole figure, ~1.5 s here, on every key press."""

    def _repaint(self, method):
        fig = self.ax.get_figure(root=True)
        canvas = fig.canvas
        if not canvas.supports_blit or fig._get_renderer() is None:
            return method()
        canvas.draw = lambda: None
        try:
            method()
        finally:
            del canvas.draw
        if not self.capturekeystrokes:
            self.cursor.set_visible(False)
        fig.draw_artist(self.ax)
        canvas.blit(self.ax.bbox)

    def _rendercursor(self):
        self._repaint(super()._rendercursor)

    def stop_typing(self):
        self._repaint(super().stop_typing)


class RefPicker:
    def __init__(self, plume_name, paths):
        self.plume_name = plume_name
        self.paths = paths           # the frames of {name}/{wavelength}/, in time order
        self.k = 0
        self.busy = False            # a frame step is running
        self.loader = ThreadPoolExecutor(max_workers=1)
        self.preloaded = {}          # frame index: Future of its map (the next one in the stepping direction)
        self.direction = +1
        m = load_map(paths[0])
        self.ax = None
        self.points = []   # [(Tx, Ty, DATE-OBS of the frame it was picked on)]
        self.artists = []

        positive = m.data[m.data > 0]
        lo, hi = np.percentile(positive, [1, 99.9])
        self.norm = ImageNormalize(vmin=lo, vmax=hi, stretch=next(iter(STRETCHES.values()))(), clip=True)

        self.fig = plt.figure(figsize=(12, 11))

        # A full redraw of the 4096x4096 image and its WCS grid takes ~1.5 s, so while a
        # slider is dragged only the slider is repainted (blitting) and the image is
        # redrawn once, when the mouse button is released.
        top = float(np.percentile(positive, 99.99))
        self.blit = self.fig.canvas.supports_blit
        self.background = None
        self.contrast_pending = False
        self.saving = False
        self.inset_mode = False      # next clicks on the image fill the inset corners
        self.corner = 0              # which corner the next click fills: 0 bottom left, 1 top right
        self.inset_ax = None
        self.inset_box = None        # rectangle on the main image
        self.box_px = {}             # its corners in image pixels, keyed (sx, sy): 0 left/bottom, 1 right/top
        self.connectors = []         # lines from the region's box to the inset
        self.drag = None             # while moving / resizing the inset, see _on_click
        self.cursor = Cursors.POINTER
        # previews while dragging (drawn by blitting only)
        self.drag_outline = Rectangle((0, 0), 0, 0, transform=self.fig.transFigure, fill=False,
                                      edgecolor=INSET_COLOR, ls="--", lw=1.5, animated=True)
        self.drag_lines = [Line2D([], [], transform=self.fig.transFigure, color=INSET_COLOR, ls="--",
                                  lw=1, animated=True) for _ in range(4)]
        for a in [self.drag_outline] + self.drag_lines:
            self.fig.add_artist(a)
        self.sliders = []
        for k, (label, val) in enumerate([("vmin", lo), ("vmax", hi)]):
            s = Slider(self.fig.add_axes([0.12, 0.08 - k * 0.04, 0.5, 0.025]), label, 0, top,
                       valinit=val, initcolor="none")
            if self.blit:
                s.drawon = False
                for a in self.slider_artists(s):
                    a.set_animated(True)
            s.on_changed(self._slider_moved)
            self.sliders.append(s)
        stretch_label = self.fig.text(0.8, 0.135, "stretch", ha="left", va="bottom")
        self.stretch_buttons = RadioButtons(self.fig.add_axes([0.8, 0.005, 0.15, 0.13], frameon=False),
                                            list(STRETCHES))
        self.stretch_buttons.on_clicked(self._stretch)
        self.save_buttons = []
        for k, fmt in enumerate(["png", "pdf"]):
            b = Button(self.fig.add_axes([0.69, 0.075 - k * 0.045, 0.09, 0.035]), f"save {fmt}",
                       color=BUTTON_COLOR, hovercolor=BUTTON_COLOR)
            b.on_clicked(lambda _, fmt=fmt: self.save(fmt))
            self.save_buttons.append(b)
        # inset row; hovercolor = color: the stock hover effect redraws the whole figure
        self.inset_button = Button(self.fig.add_axes([0.04, 0.125, 0.09, 0.035]), "inset",
                                   color=BUTTON_COLOR, hovercolor=BUTTON_COLOR)
        self.inset_button.on_clicked(self._toggle_inset_mode)
        self.corner_boxes = []
        for label, x0 in [("bottom left", 0.25), ("top right", 0.5)]:
            tb = QuickTextBox(self.fig.add_axes([x0, 0.125, 0.15, 0.035]), label,
                              color="white", hovercolor="white")
            tb.on_submit(lambda _: self.make_inset())
            self.corner_boxes.append(tb)
        self.remove_button = Button(self.fig.add_axes([0.69, 0.125, 0.09, 0.035]), "remove inset",
                                    color=BUTTON_COLOR, hovercolor=BUTTON_COLOR)
        self.remove_button.on_clicked(lambda _: self.remove_inset(clear_corners=True))
        # everything hidden when saving
        self.widgets = ([s.ax for s in self.sliders] + [self.stretch_buttons.ax, stretch_label]
                        + [b.ax for b in self.save_buttons]
                        + [self.inset_button.ax, self.remove_button.ax] + [tb.ax for tb in self.corner_boxes])

        # spacecraft views (right panel)
        self.traj = Trajectories()
        self.view_artists = []
        title = self.fig.text(0.8, 0.92, "Spacecraft view\n(JPL Horizons at DATE-OBS)", va="top", fontsize=10)
        names = list(VIEWS)
        self.view_buttons = CheckButtons(self.fig.add_axes([0.79, 0.715, 0.2, 0.17], frameon=False), names,
                                         [n in VIEWS_SHOWN for n in names],
                                         check_props={"facecolor": list(VIEWS.values()), "linewidth": 2})
        self.view_buttons.on_clicked(lambda _: self.draw_views())
        self.view_info = self.fig.text(0.8, 0.705, "", va="top", fontsize=8, family="monospace")
        self.v_range, self.margin = V_RANGE, MARGIN
        self.connect_toggle = CheckButtons(self.fig.add_axes([0.79, 0.41, 0.2, 0.035], frameon=False),
                                           ["connectivity region"], [True])
        self.connect_toggle.on_clicked(lambda _: self.draw_views())
        self.speed_box = QuickTextBox(self.fig.add_axes([0.845, 0.38, 0.075, 0.026]), "km/s ",
                                      initial=f"{V_RANGE[0]}-{V_RANGE[1]}", color="white", hovercolor="white")
        self.speed_box.on_submit(self._connect_typed)
        self.margin_box = QuickTextBox(self.fig.add_axes([0.955, 0.38, 0.035, 0.026]), "± ",
                                       initial=f"{MARGIN}", color="white", hovercolor="white")
        self.margin_box.on_submit(self._connect_typed)
        self.region_style = REGION_STYLES[0]
        self.style_button = Button(self.fig.add_axes([0.8, 0.348, 0.19, 0.028]),
                                   f"style: {self.region_style}", color=BUTTON_COLOR, hovercolor=BUTTON_COLOR)
        self.style_button.on_clicked(self._next_style)
        key = self.fig.text(0.8, 0.268, "x  sub-spacecraft point\n   (centre of its disk)\n"
                                        "o  same, on the far side\n"
                                        "—  its limb; dotted: hidden\n"
                                        "▭  connectivity region:\n   Parker spiral, no PFSS",
                            va="bottom", fontsize=8, family="monospace", color="0.3")
        self.widgets += [self.view_buttons.ax, title, self.view_info, key, self.connect_toggle.ax,
                         self.speed_box.ax, self.margin_box.ax, self.style_button.ax]

        # frames (only with several)
        self.frame_buttons = []
        if len(paths) > 1:
            self.frame_label = self.fig.text(0.8, 0.235, "", fontsize=9, va="bottom")
            for label, x0, step in [("< prev", 0.8, -1), ("next >", 0.89, +1)]:
                b = Button(self.fig.add_axes([x0, 0.19, 0.08, 0.035]), label,
                           color=BUTTON_COLOR, hovercolor=BUTTON_COLOR)
                b.on_clicked(lambda _, step=step: self.step(step))
                self.frame_buttons.append(b)
            self.widgets += [self.frame_label] + [b.ax for b in self.frame_buttons]

        self.show_frame(0, m, redraw=False)
        self.fig.canvas.mpl_connect("draw_event", self._on_draw)
        self.fig.canvas.mpl_connect("button_release_event", self._on_release)
        self.fig.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.fig.canvas.mpl_connect("button_press_event", self._on_click)
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self.fig.canvas.mpl_connect("close_event", self._on_close)

    # ---- frames ----
    def show_frame(self, k, m=None, redraw=True):
        """Show frame k: a new image axes with that frame's coordinates, and everything drawn
        on it again (grid, labels, points, inset, spacecraft views)."""
        self.k = k
        path = self.paths[k]
        if m is None:
            future = self.preloaded.pop(k, None)
            m = future.result() if future is not None else load_map(path)
        self.m = m
        self.save_stem = path.with_suffix("")   # for "save png" / "save pdf"
        limits = inset_pos = None
        if self.ax is not None:
            limits = self.ax.get_xlim(), self.ax.get_ylim()
            if self.inset_ax is not None:
                inset_pos = self.inset_ax.get_position()
            self.remove_inset()
            self.ax.remove()
        self.artists, self.view_artists = [], []

        self.ax = self.fig.add_axes(AX_POS, projection=m)
        self.image = m.plot(axes=self.ax, norm=self.norm)
        # shrink the 4096x4096 data to the screen first, then apply the colour scale: by default
        # matplotlib colours all 16.8 M pixels first (~1 GB of temporary memory, 1.2 s per redraw)
        self.image.set_interpolation_stage("data")
        self.image.norm = self.norm  # sunpy plots with a copy of the norm; share ours (sliders, stretch)
        m.draw_limb(axes=self.ax, color="white", lw=0.8)
        grid = m.draw_grid(axes=self.ax, color="white", alpha=0.3, lw=0.5)
        # Stonyhurst longitude ticks on the top only, latitude on the right only (sunpy allows
        # both on both, so near the poles longitudes could appear on the right)
        for coord, side in ((grid[0], "t"), (grid[1], "r")):
            coord.set_ticks_position(side)
            coord.set_ticklabel_position(side)
            coord.set_axislabel_position(side)
        self.label_grid()
        # above the Stonyhurst longitude ticks of the top axis
        self.ax.set_title(f"{self.plume_name}   AIA {int(m.wavelength.value)} Å   {m.date.isot[:19]}", pad=32)
        obs = m.observer_coordinate
        self.observer = (m.dsun.to_value(u.m), m.rsun_meters.to_value(u.m),
                         obs.lat.to_value(u.rad), obs.lon.to_value(u.rad))
        self.ax.format_coord = lambda x, y: self.format_coord(self.m, x, y)
        if limits is not None:   # keep the zoom (frames have nearly the same pointing)
            self.ax.set_xlim(*limits[0])
            self.ax.set_ylim(*limits[1])
        self.ax.callbacks.connect("xlim_changed", lambda _: self.update_connectors())
        self.ax.callbacks.connect("ylim_changed", lambda _: self.update_connectors())

        self.draw(redraw=False)
        if inset_pos is not None:
            self.ax.apply_aspect()
            self.make_inset(position=inset_pos)
        self.draw_views(redraw=False)
        n = len(self.paths)
        self.fig.canvas.manager.set_window_title(
            f"Reference file: {self.plume_name}" + (f"   frame {k + 1}/{n}" if n > 1 else ""))
        if n > 1:
            self.frame_label.set_text(f"frame {k + 1} / {n}   {m.date.isot[:19]}")
            # reading a frame takes ~2 s: load the next one (in the stepping direction) while this
            # one is shown; only one, to keep memory low
            keep = {k + self.direction} & set(range(n))
            for j in list(self.preloaded):
                if j not in keep:
                    self.preloaded.pop(j).cancel()
            for j in keep - set(self.preloaded):
                self.preloaded[j] = self.loader.submit(load_map, self.paths[j])
        if redraw:
            self.fig.canvas.draw_idle()

    def step(self, d):
        if self.busy:   # clicks / keys that arrive during a step are dropped (see below)
            return
        k = min(max(self.k + d, 0), len(self.paths) - 1)
        self.direction = +1 if d > 0 else -1
        if k == self.k:
            return
        self.busy = True
        manager = self.fig.canvas.manager
        try:
            manager.set_window_title(f"Reference file: {self.plume_name}   loading frame {k + 1}/{len(self.paths)} ...")
            self.show_frame(k, redraw=False)
            # the previous frame's map and axes sit in reference cycles, which Python frees only
            # now and then: without this, memory grows by ~70 MB per step before it does
            gc.collect()
            self.fig.canvas.draw()
            # handle the events that queued up meanwhile now, while busy: otherwise every extra
            # click on next (or a held n key) would replay a full step afterwards
            self.fig.canvas.flush_events()
        finally:
            self.busy = False

    # ---- spacecraft views ----
    def _connect_typed(self, _):
        try:
            v_range, margin = parse_speeds(self.speed_box.text), float(self.margin_box.text)
        except ValueError:
            print(f"Connectivity: could not read {self.speed_box.text!r} km/s, ± {self.margin_box.text!r} deg "
                  "(e.g. 300-800 and 10)")
            return
        if (v_range, margin) != (self.v_range, self.margin):
            self.v_range, self.margin = v_range, margin
            self.draw_views()

    def _next_style(self, _):
        self.region_style = REGION_STYLES[(REGION_STYLES.index(self.region_style) + 1) % len(REGION_STYLES)]
        self.style_button.label.set_text(f"style: {self.region_style}")
        self.draw_views()

    def dim_outside(self, regions):
        """Darken the image outside the union of the regions (polygons in image pixels)."""
        ny, nx = self.m.data.shape
        step = 4   # mask resolution: 1024 x 1024 for a full-disk image
        gx, gy = np.meshgrid(np.arange(0, nx, step) + step / 2 - 0.5, np.arange(0, ny, step) + step / 2 - 0.5)
        points = np.c_[gx.ravel(), gy.ravel()]
        inside = np.zeros(len(points), bool)
        for xy in regions:
            inside |= MplPath(xy).contains_points(points)
        alpha = np.where(inside.reshape(gx.shape), 0.0, DIM)
        rgba = np.zeros(gx.shape + (4,))
        rgba[..., 3] = alpha
        limits = self.ax.get_xlim(), self.ax.get_ylim()   # imshow could reset the zoom
        # above the AIA image (zorder 0), below the outlines, labels and grid
        self.view_artists.append(self.ax.imshow(rgba, origin="lower", extent=(-0.5, nx - 0.5, -0.5, ny - 0.5),
                                                interpolation="bilinear", zorder=0.5))
        self.ax.set_xlim(*limits[0])
        self.ax.set_ylim(*limits[1])

    def visible(self, lon, lat):
        """True where Stonyhurst lon, lat (deg, arrays) on the surface face the AIA observer."""
        obs = self.m.observer_coordinate
        o_lon, o_lat = obs.lon.to_value(u.rad), obs.lat.to_value(u.rad)
        lon, lat = np.radians(lon), np.radians(lat)
        cos_sep = np.sin(lat) * np.sin(o_lat) + np.cos(lat) * np.cos(o_lat) * np.cos(lon - o_lon)
        return cos_sep > (self.m.rsun_meters / self.m.dsun).decompose().value

    def surface_to_pixels(self, lon, lat):
        """Stonyhurst lon, lat (deg, arrays) on the surface -> image pixels; points on the far side
        are moved to the limb, in their direction."""
        m = self.m
        hpc = SkyCoord(lon * u.deg, lat * u.deg, m.rsun_meters,
                       frame=HeliographicStonyhurst(obstime=m.date)).transform_to(m.coordinate_frame)
        tx, ty = hpc.Tx.to_value(u.arcsec), hpc.Ty.to_value(u.arcsec)
        hidden = ~self.visible(lon, lat)
        if hidden.any():
            r_limb = m.rsun_obs.to_value(u.arcsec)
            angle = np.arctan2(ty[hidden], tx[hidden])
            tx[hidden], ty[hidden] = r_limb * np.cos(angle), r_limb * np.sin(angle)
        px = m.world_to_pixel(SkyCoord(tx * u.arcsec, ty * u.arcsec, frame=m.coordinate_frame))
        return px.x.value, px.y.value

    def draw_connectivity(self, c, name, color):
        """The connectivity region of a spacecraft at c (see the docstring); returns the
        Stonyhurst longitudes reached for the fast and slow speed, and the region outline in
        image pixels (None if all on the far side)."""
        lon0, lat0 = c.lon.to_value(u.deg), c.lat.to_value(u.deg)
        fast, slow = (lon0 + spiral_shift(c.radius, R_SS, v) for v in self.v_range[::-1])
        lon_a, lon_b = fast - self.margin, slow + self.margin
        lat_a, lat_b = max(lat0 - self.margin, -89.9), min(lat0 + self.margin, 89.9)
        n = 60   # points per edge, so the edges follow the sphere
        lon = np.r_[np.linspace(lon_a, lon_b, n), np.full(n, lon_b), np.linspace(lon_b, lon_a, n), np.full(n, lon_a)]
        lat = np.r_[np.full(n, lat_a), np.linspace(lat_a, lat_b, n), np.full(n, lat_b), np.linspace(lat_b, lat_a, n)]
        if not self.visible(lon, lat).any() and not self.visible(np.array([fast, slow]), np.array([lat0] * 2)).any():
            return fast, slow, None   # all on the far side
        x, y = self.surface_to_pixels(lon, lat)
        region = np.c_[x, y]
        if self.region_style == "shaded":
            patch = Polygon(region, closed=True, facecolor=to_rgba(color, 0.18), edgecolor=color, lw=1)
        else:   # outline: coloured, with a dark edge so it shows on bright and dark parts alike
            patch = Polygon(region, closed=True, fill=False, edgecolor=color, lw=2, path_effects=OUTLINE)
        self.ax.add_patch(patch)
        self.view_artists.append(patch)
        # name at the middle of the region's north edge
        top = 2 * n + n // 2
        self.view_artists.append(self.ax.annotate(name, (x[top], y[top]), xytext=(0, 4), textcoords="offset points",
                                                  ha="center", va="bottom", color=color, fontsize=9,
                                                  path_effects=OUTLINE))
        # the speed range itself: a line along the spacecraft's latitude, ends labelled in km/s
        line_lon = np.linspace(fast, slow, 40)
        x, y = self.surface_to_pixels(line_lon, np.full_like(line_lon, lat0))
        self.view_artists += self.ax.plot(x, y, color=color, lw=2.5, path_effects=OUTLINE)
        for (xe, ye), v, ha in (((x[0], y[0]), self.v_range[1], "right"), ((x[-1], y[-1]), self.v_range[0], "left")):
            self.view_artists.append(self.ax.annotate(f"{v:g}", (xe, ye), xytext=(0, -9), textcoords="offset points",
                                                      ha=ha, color=color, fontsize=8, path_effects=OUTLINE))
        return fast, slow, region

    def draw_views(self, redraw=True):
        for a in self.view_artists:
            a.remove()
        self.view_artists = []
        m, t = self.m, self.m.date
        shown = [n for n, on in zip(VIEWS, self.view_buttons.get_status()) if on]
        self.traj.fetch(t, [n for n in shown if n != L1], status=print)
        connect = self.connect_toggle.get_status()[0]
        carrington = carrington_L0(t).to_value(u.deg)   # Carrington longitude = Stonyhurst + L0
        info, regions = [], []
        for name in shown:
            if name == L1:
                earth = get_body_heliographic_stonyhurst("earth", t)
                c = HeliographicStonyhurst(earth.lon, earth.lat, 0.99 * earth.radius, obstime=t)
            else:
                c = self.traj.position(name, t)
            if c is None:
                info.append(f"{name}: not available")
                continue
            color = VIEWS[name]
            lines = [name]
            if connect:
                fast, slow, region = self.draw_connectivity(c, name, color)
                if region is not None:
                    regions.append(region)
                lines.append(f"  connects lon {fast:+.0f}..{slow:+.0f}°")
                lines.append(f"   (Carr. {(fast + carrington) % 360:.0f}..{(slow + carrington) % 360:.0f}°) ±{self.margin:g}°")
            if name == L1:
                info.append("\n".join(lines))
                continue
            obs = m.observer_coordinate
            o_lon, o_lat = obs.lon.to_value(u.rad), obs.lat.to_value(u.rad)
            lon, lat = c.lon.to_value(u.rad), c.lat.to_value(u.rad)
            # angle at the Sun between the spacecraft and the AIA observer
            cos_sep = np.sin(lat) * np.sin(o_lat) + np.cos(lat) * np.cos(o_lat) * np.cos(lon - o_lon)
            sep = np.degrees(np.arccos(np.clip(cos_sep, -1, 1)))
            on_disk = cos_sep > (m.rsun_meters / m.dsun).decompose().value
            subpoint = SkyCoord(c.lon, c.lat, m.rsun_meters, frame=HeliographicStonyhurst(obstime=t))
            if on_disk:
                px = m.world_to_pixel(subpoint)
                marker, label = "x", name
            else:   # at the limb, in the direction of the hidden sub-point
                hpc = subpoint.transform_to(m.coordinate_frame)
                angle = np.arctan2(hpc.Ty.to_value(u.arcsec), hpc.Tx.to_value(u.arcsec))
                r = m.rsun_obs.to_value(u.arcsec)
                px = m.world_to_pixel(SkyCoord(r * np.cos(angle) * u.arcsec, r * np.sin(angle) * u.arcsec,
                                               frame=m.coordinate_frame))
                marker, label = "o", f"{name} (far side)"
            x, y = float(px.x.value), float(px.y.value)
            self.view_artists += self.ax.plot(x, y, marker, color=color, ms=10, mew=2,
                                              mfc="none", path_effects=OUTLINE)
            self.view_artists.append(self.ax.annotate(label, (x, y), xytext=(7, 5), textcoords="offset points",
                                                      color=color, fontsize=9, path_effects=OUTLINE))
            visible, hidden = limb(self.ax, SkyCoord(c), rsun=m.rsun_meters, color=color, lw=1.3)
            self.view_artists += [a for a in (visible, hidden) if a is not None]
            lines[1:1] = [f"  r {c.radius.to_value(u.R_sun):.1f} Rsun, {sep:.1f}° from SDO",
                          f"  lon {c.lon.wrap_at(180 * u.deg).to_value(u.deg):+.1f} lat {c.lat.to_value(u.deg):+.1f}"
                          f" ({'on disk' if on_disk else 'far side'})"]
            info.append("\n".join(lines))
        if connect and regions and self.region_style == "dim outside":
            self.dim_outside(regions)
        self.view_info.set_text("\n".join(info))
        if redraw:
            self.fig.canvas.draw_idle()

    def label_grid(self, spacing=15):
        """Label the heliographic grid on the image (the axes ticks only show where grid lines
        reach the frame, i.e. never on a full-disk view): latitudes just outside both limbs,
        longitudes along the equator."""
        obs = self.m.observer_coordinate
        L0, B0 = obs.lon.to_value(u.deg), obs.lat.to_value(u.deg)
        frame = HeliographicStonyhurst(obstime=self.m.date)
        kw = dict(color="white", fontsize=8, textcoords="offset points", path_effects=OUTLINE)
        labels = []   # (lon, lat, text, offset, ha, va); converted to pixels in one go below
        for lat in range(-90 + spacing, 90, spacing):
            # the limb at this latitude: where the point is 90 deg from the observer
            dlon = np.degrees(np.arccos(-np.tan(np.radians(lat)) * np.tan(np.radians(B0))))
            for side in (+1, -1):   # west, east limb
                labels.append((L0 + side * dlon, lat, f"{lat}°", (side * 4, 0),
                               "left" if side > 0 else "right", "center"))
        for lon in range(-180, 180, spacing):
            if np.cos(np.radians(lon - L0)) > np.cos(np.radians(80)):   # on the visible side
                labels.append((lon, 0, f"{lon}°", (2, -2), "left", "top"))
        lons, lats = (np.array([lab[k] for lab in labels], dtype=float) for k in (0, 1))
        px = self.m.world_to_pixel(SkyCoord(lons * u.deg, lats * u.deg, frame=frame))
        for x, y, (_, _, text, offset, ha, va) in zip(px.x.value, px.y.value, labels):
            self.ax.annotate(text, (x, y), xytext=offset, ha=ha, va=va, **kw)

    def hgs(self, tx, ty):
        """Helioprojective (arcsec) -> Stonyhurst lon, lat (deg) on the solar surface, or None
        off disk. Closed form (Thompson 2006, A&A 449, 791): the sunpy transform takes ~40 ms,
        too slow for every mouse move."""
        D, R, B0, L0 = self.observer
        tx, ty = np.radians(tx / 3600), np.radians(ty / 3600)
        q = D * np.cos(ty) * np.cos(tx)
        disc = q ** 2 - D ** 2 + R ** 2
        if disc < 0:
            return None
        d = q - np.sqrt(disc)   # distance observer -> surface
        x, y, z = d * np.cos(ty) * np.sin(tx), d * np.sin(ty), D - d * np.cos(ty) * np.cos(tx)
        lat = np.arcsin((y * np.cos(B0) + z * np.sin(B0)) / R)
        lon = L0 + np.arctan2(x, z * np.cos(B0) - y * np.sin(B0))
        return (np.degrees(lon) + 180) % 360 - 180, np.degrees(lat)

    def format_coord(self, m, x, y):
        tx, ty = (3600 * ((v + 180) % 360 - 180) for v in m.wcs.pixel_to_world_values(x, y))
        hgs = self.hgs(tx, ty)
        hg = f"lon {hgs[0]:.2f}°  lat {hgs[1]:.2f}°" if hgs else "off disk"
        return f"x {tx:.1f}\"  y {ty:.1f}\"   |   HGS {hg}"

    def add_point(self, x_px, y_px):
        c = self.m.pixel_to_world(x_px * u.pix, y_px * u.pix)
        tx, ty = c.Tx.to_value(u.arcsec), c.Ty.to_value(u.arcsec)
        self.points.append((tx, ty, self.m.date.isot[:19]))
        print(f"point {len(self.points)}: x={tx:.1f}\"  y={ty:.1f}\"  ({self.m.date.isot[:19]})")
        self.draw()

    def draw(self, redraw=True):
        for a in self.artists:
            a.remove()
        self.artists = []
        for k, (tx, ty, _) in enumerate(self.points, 1):
            px = self.m.world_to_pixel(SkyCoord(tx * u.arcsec, ty * u.arcsec, frame=self.m.coordinate_frame))
            x, y = float(px.x.value), float(px.y.value)
            last = k == len(self.points)
            color = "red" if last else "cyan"
            self.artists += self.ax.plot(x, y, "+", color=color, ms=14, mew=1.5)
            self.artists.append(self.ax.text(x, y, f"  {k}: ({tx:.1f}\", {ty:.1f}\")", color=color,
                                             fontsize=10, va="bottom"))
        if redraw:
            self.fig.canvas.draw_idle()

    @staticmethod
    def slider_artists(s):
        return [s.poly, s._handle, s.valtext]

    def draw_sliders(self):
        for s in self.sliders:
            for a in self.slider_artists(s):
                self.fig.draw_artist(a)

    def repaint_over_background(self, *artists):
        # restore the last full redraw, then repaint what may have changed since (cheap)
        self.fig.canvas.restore_region(self.background)
        for a in [self.inset_button.ax] + [tb.ax for tb in self.corner_boxes]:
            self.fig.draw_artist(a)
        self.draw_sliders()
        for a in artists:
            self.fig.draw_artist(a)
        self.fig.canvas.blit(self.fig.bbox)

    def blit_axes(self, ax):
        if self.blit and self.background is not None:
            self.fig.draw_artist(ax)
            self.fig.canvas.blit(ax.bbox)
        else:
            self.fig.canvas.draw_idle()

    # ---- inset ----
    def _toggle_inset_mode(self, _=None):
        self.inset_mode = not self.inset_mode
        self.corner = 0
        b = self.inset_button
        b.color = b.hovercolor = ACTIVE_COLOR if self.inset_mode else BUTTON_COLOR
        b.ax.set_facecolor(b.color)
        b.label.set_text("click 2 corners" if self.inset_mode else "inset")
        self.blit_axes(b.ax)

    def set_corner(self, k, text, notify=True):
        tb = self.corner_boxes[k]
        tb.eventson = notify
        tb.set_val(text)
        tb.eventson = True

    def corner_coord(self, k):
        text = self.corner_boxes[k].text.replace(",", " ").split()
        try:
            x, y = (float(v) for v in text)
        except ValueError:
            return None
        return x, y

    def make_inset(self, position=None):
        corners = [self.corner_coord(k) for k in range(2)]
        if None in corners:
            return
        (x0, y0), (x1, y1) = corners
        x0, x1, y0, y1 = min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1)
        if x0 == x1 or y0 == y1:
            print("Inset: the two corners must differ in x and in y")
            return
        if self.inset_ax is not None and (x0, x1, y0, y1) == self.inset_region:
            return
        bl = SkyCoord(x0 * u.arcsec, y0 * u.arcsec, frame=self.m.coordinate_frame)
        tr = SkyCoord(x1 * u.arcsec, y1 * u.arcsec, frame=self.m.coordinate_frame)
        sub = self.m.submap(bl, top_right=tr)
        if min(sub.data.shape) < 2:
            print("Inset: region too small or outside the image")
            return

        old = self.inset_ax.get_position() if self.inset_ax is not None else position
        self.remove_inset()
        # 30% of the image width, same shape as the region; a new region keeps the
        # inset's place (top left corner) and width
        pos = self.ax.get_position()
        fw, fh = self.fig.get_size_inches()
        ny, nx = sub.data.shape
        w = old.width if old is not None else 0.3 * pos.width
        h = w * fw / fh * ny / nx
        left, top = (old.x0, old.y1) if old is not None else (pos.x1 - w - 0.01, pos.y1 - 0.01)
        self.inset_ax = self.fig.add_axes([left, top - h, w, h], projection=sub)
        self.inset_image = sub.plot(axes=self.inset_ax, title=False)
        self.inset_image.set_interpolation_stage("data")
        self.inset_image.norm = self.norm   # follow the sliders and the stretch
        self.inset_ax.grid(False)   # sunpy adds a dotted coordinate grid
        for c in (self.inset_ax.coords[0], self.inset_ax.coords[1]):
            c.set_ticks_visible(False)
            c.set_ticklabel_visible(False)
            c.set_axislabel("")
        self.inset_ax.coords.frame.set_color(INSET_COLOR)
        self.inset_ax.coords.frame.set_linewidth(1.5)
        self.inset_ax.format_coord = lambda x, y, sub=sub: self.format_coord(sub, x, y)
        self.inset_box = self.m.draw_quadrangle(bl, top_right=tr, axes=self.ax,
                                                edgecolor=INSET_COLOR, lw=1.2)
        for sx, x in enumerate((x0, x1)):
            for sy, y in enumerate((y0, y1)):
                px = self.m.world_to_pixel(SkyCoord(x * u.arcsec, y * u.arcsec, frame=self.m.coordinate_frame))
                self.box_px[sx, sy] = (float(px.x.value), float(px.y.value))
        self.update_connectors()
        self.inset_region = (x0, x1, y0, y1)
        print(f"Inset: x {x0:.1f}\" to {x1:.1f}\", y {y0:.1f}\" to {y1:.1f}\"")
        if self.inset_mode:
            self._toggle_inset_mode()
        self.fig.canvas.draw_idle()

    def remove_inset(self, clear_corners=False):
        if self.inset_ax is not None:
            self.inset_ax.remove()
            self.inset_box.remove()
            self.inset_ax = self.inset_box = None
        self.update_connectors()
        if clear_corners:
            for k in range(2):
                self.set_corner(k, "", notify=False)
            self.fig.canvas.draw_idle()

    def _corner_click(self, x_px, y_px):
        c = self.m.pixel_to_world(x_px * u.pix, y_px * u.pix)
        text = f"{c.Tx.to_value(u.arcsec):.1f}, {c.Ty.to_value(u.arcsec):.1f}"
        if self.corner == 0:
            self.set_corner(1, "", notify=False)
            self.set_corner(0, text, notify=False)
            self.corner = 1
        else:
            self.set_corner(1, text)   # submit -> make_inset

    def connector_pairs(self, bounds):
        """Corner pairs ((sx, sy) of the box, (sx, sy) of the inset) to join with lines: the
        edges of the convex hull of both rectangles that go from one to the other, i.e.
        the outer lines, which cross neither rectangle. bounds: inset in figure fraction."""
        to_fig = self.fig.transFigure.inverted()
        x0, y0, w, h = bounds
        pts = [(*to_fig.transform(self.ax.transData.transform(xy)), ("box", k)) for k, xy in self.box_px.items()]
        pts += [(x0 + sx * w, y0 + sy * h, ("inset", (sx, sy))) for sx in (0, 1) for sy in (0, 1)]
        hull = convex_hull(pts)
        pairs = []
        for a, b in zip(hull, hull[1:] + hull[:1]):
            if a[2][0] != b[2][0]:
                box, inset = (a, b) if a[2][0] == "box" else (b, a)
                pairs.append(((box[0], box[1]), box[2][1], (inset[0], inset[1]), inset[2][1]))
        return pairs   # [(box xy in figure fraction, box key, inset xy in figure fraction, inset key)]

    def update_connectors(self):
        for c in self.connectors:
            c.remove()
        self.connectors = []
        if self.inset_ax is None:
            return
        for _, kb, _, ki in self.connector_pairs(self.inset_ax.get_position().bounds):
            c = ConnectionPatch(xyA=self.box_px[kb], coordsA="data", axesA=self.ax,
                                xyB=ki, coordsB="axes fraction", axesB=self.inset_ax,
                                color=INSET_COLOR, lw=1)
            self.fig.add_artist(c)
            self.connectors.append(c)

    def inset_hit(self, event):
        """None, "move" (inside the inset) or the (sx, sy) corner under the mouse."""
        if self.inset_ax is None:
            return None
        b = self.inset_ax.bbox
        for sx, cx in enumerate((b.x0, b.x1)):
            for sy, cy in enumerate((b.y0, b.y1)):
                if abs(event.x - cx) <= HANDLE_PX and abs(event.y - cy) <= HANDLE_PX:
                    return sx, sy
        return "move" if b.contains(event.x, event.y) else None

    def set_cursor(self, cursor):
        if cursor != self.cursor:
            self.cursor = cursor
            self.fig.canvas.set_cursor(cursor)

    def _on_motion(self, event):
        if self.drag is None:
            toolbar = getattr(self.fig.canvas.manager, "toolbar", None)
            if not (toolbar is not None and getattr(toolbar, "mode", "")):
                hit = self.inset_hit(event)
                self.set_cursor(Cursors.POINTER if hit is None else
                                Cursors.MOVE if hit == "move" else Cursors.SELECT_REGION)
            return
        d = self.drag
        fw, fh = self.fig.bbox.width, self.fig.bbox.height
        b = d["bbox"]   # inset at the press, display pixels
        if d["corner"] is None:   # move
            x0, y0, w, h = b.x0 + event.x - d["x"], b.y0 + event.y - d["y"], b.width, b.height
        else:                     # resize, keeping the shape and the opposite corner
            sx, sy = d["corner"]
            ax_, ay_ = (b.x0, b.x1)[1 - sx], (b.y0, b.y1)[1 - sy]
            aspect = b.height / b.width
            w = max(abs(event.x - ax_), abs(event.y - ay_) / aspect, MIN_INSET_PX)
            h = w * aspect
            x0 = ax_ if sx else ax_ - w
            y0 = ay_ if sy else ay_ - h
        bounds = (x0 / fw, y0 / fh, w / fw, h / fh)
        self.drag_outline.set_bounds(*bounds)
        pairs = self.connector_pairs(bounds)
        for k, line in enumerate(self.drag_lines):
            if k < len(pairs):
                (bx, by), _, (ix, iy), _ = pairs[k]
                line.set_data([bx, ix], [by, iy])
            else:
                line.set_data([], [])
        if self.blit and self.background is not None:
            self.repaint_over_background(self.drag_outline, *self.drag_lines)

    def save(self, fmt):
        out = self.save_stem.with_suffix(f".{fmt}")
        self.saving = True
        for w in self.widgets:
            w.set_visible(False)
        try:
            self.fig.savefig(out, dpi=150, bbox_inches="tight")
        finally:
            for w in self.widgets:
                w.set_visible(True)
            self.saving = False
            self.fig.canvas.draw_idle()
        print(f"Saved {out}")

    def _on_draw(self, _):
        # after each full redraw on screen: keep the figure without the sliders' moving parts
        # (not when saving a file: savefig also sends draw events, and may use a pdf canvas)
        if self.blit and not self.saving and self.fig.canvas.supports_blit:
            self.background = self.fig.canvas.copy_from_bbox(self.fig.bbox)
            self.draw_sliders()

    def _slider_moved(self, _):
        self.contrast_pending = True
        if self.blit and self.background is not None:
            self.repaint_over_background()

    def images_changed(self):
        self.image.changed()
        if self.inset_ax is not None:
            self.inset_image.changed()
        self.fig.canvas.draw_idle()

    def _stretch(self, name):
        self.norm.stretch = STRETCHES[name]()
        self.norm.inverse_stretch = self.norm.stretch.inverse
        self.images_changed()

    def _on_release(self, _):
        if self.drag is not None:
            self.drag = None
            self.inset_ax.set_position(self.drag_outline.get_bbox().bounds)
            self.update_connectors()
            self.fig.canvas.draw_idle()
        if not self.contrast_pending:
            return
        self.contrast_pending = False
        vmin, vmax = (s.val for s in self.sliders)
        if vmin < vmax:
            self.norm.vmin, self.norm.vmax = vmin, vmax
        self.images_changed()

    def _on_click(self, event):
        toolbar = getattr(self.fig.canvas.manager, "toolbar", None)
        if toolbar is not None and getattr(toolbar, "mode", ""):
            return
        if event.button != 1:
            return
        hit = self.inset_hit(event)
        if hit is not None:
            self.drag = {"x": event.x, "y": event.y, "bbox": self.inset_ax.bbox.frozen(),
                         "corner": None if hit == "move" else hit}
            self.drag_outline.set_bounds(*self.inset_ax.get_position().bounds)
        elif event.inaxes is self.ax:
            if self.inset_mode:
                self._corner_click(event.xdata, event.ydata)
            else:
                self.add_point(event.xdata, event.ydata)

    def _on_key(self, event):
        if any(tb.capturekeystrokes for tb in self.corner_boxes + [self.speed_box, self.margin_box]):
            return
        if event.key == "u" and self.points:
            self.points.pop()
            self.draw()
        elif event.key in ("n", "b") and len(self.paths) > 1:
            self.step(+1 if event.key == "n" else -1)

    def _on_close(self, event):
        if self.points:
            tx, ty, date_obs = self.points[-1]
            save_coord(self.plume_name, date_obs, tx, ty)
        else:
            print("No point picked, nothing saved")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--name", required=True, help="e.g. plume_002: shows {name}/{wavelength}/")
    parser.add_argument("--wavelength", type=int, default=171, choices=WAVELENGTHS,
                        help="AIA channel in Å (default 171)")
    args = parser.parse_args()

    folder = data_dir(args.name, args.wavelength)
    paths = sorted(folder.glob("*.fits"), key=frame_time) if folder.is_dir() else []
    if not paths:
        old = OUT_DIR / (f"{args.name}.fits" if args.wavelength == 171 else f"{args.name}_{args.wavelength}.fits")
        hint = (f"\n{old.name} is in the old layout; move it with:  mkdir -p {folder.relative_to(OUT_DIR)} "
                f"&& mv {old.name} {folder.relative_to(OUT_DIR)}/" if old.is_file() else "")
        raise SystemExit(f"No .fits in {folder} -- run download_ref_file.py first{hint}")
    print(f"{len(paths)} frame(s) in {folder}")
    picker = RefPicker(args.name, paths)  # keep a reference: matplotlib holds callbacks only weakly
    plt.show()


if __name__ == "__main__":
    main()
