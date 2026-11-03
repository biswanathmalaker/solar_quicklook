#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
Positions of the inner planets and of solar spacecraft (PSP, Solar Orbiter,
STEREO-A, SDO, SOHO, BepiColombo) on a polar plot in Stonyhurst heliographic
coordinates (Sun at the centre, Earth always at longitude 0), at a date you
choose and step through.

  date box       type a date (UTC, e.g. 2022-12-11 13:16) and press Enter
  << / >>        step back / forward by the selected step
  step           1h, 4h, 12h, 1d, 3d, 1wk, 1mo (= 30 days)
  now            go to the current time
  table          tick / untick a row to show / hide that object (default: Sun,
                 Earth, Solar Orbiter, PSP)
  Parker spirals tick box + speed range (default 300-800 km/s): the field lines from
                 each shown spacecraft (and Earth) down to the source surface
                 (2.5 Rsun) for the slowest and fastest speed, the band between
                 them shaded: where on the Sun it can be magnetically connected
  save png/pdf   save a figure of the ticked objects only (plot scaled to the
                 farthest of them, their table rows, legend) next to this script
                 as positions_YYYYMMDD_HHMM.png / .pdf

The table on the right gives each body's distance from the Sun (AU and solar
radii) and its Stonyhurst longitude / latitude (the longitude is the angle from
Earth). Spacecraft are downloaded only when ticked.

Planets come from the astropy ephemeris (no network). Spacecraft come from JPL
Horizons: their trajectories are downloaded in 120-day blocks (about 4 s per
spacecraft; Horizons refuses parallel requests) and interpolated, so steps
within a block are instant. Every block downloaded is added to
cache/horizons_positions.json in this folder (heliocentric inertial x, y, z in
AU every 2 h), which is read first, so a period is downloaded only once
(predictions, i.e. blocks reaching past the date they were downloaded, are
downloaded again after 30 days). A spacecraft is left out outside its mission
(or when Horizons has no data for that date).

Run:  ./plot_positions.py
      ./plot_positions.py --date "2022-12-11 13:16" --step 1d
"""
import argparse
import json
import math
import time
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import CartesianRepresentation
from astropy.time import Time
from matplotlib.figure import Figure
from matplotlib.widgets import Button, CheckButtons, RadioButtons, TextBox
from sunpy.coordinates import (HeliocentricInertial, HeliographicStonyhurst,
                               get_body_heliographic_stonyhurst, get_horizons_coord)
from sunpy.sun.constants import sidereal_rotation_rate
from sunpy.time import parse_time

STEPS = {"1h": 1 * u.h, "4h": 4 * u.h, "12h": 12 * u.h, "1d": 1 * u.day, "3d": 3 * u.day,
         "1wk": 7 * u.day, "1mo": 30 * u.day}
PLANETS = {   # name: (ephemeris body, colour, marker)
    "Mercury": ("mercury", "0.5", "o"),
    "Venus": ("venus", "goldenrod", "o"),
    "Earth": ("earth", "tab:blue", "o"),
    "Mars": ("mars", "tab:red", "o"),
}
SPACECRAFT = {   # name: (JPL Horizons name, first date with data, colour, marker)
    "PSP": ("Parker Solar Probe", "2018-08-13", "purple", "D"),
    "Solar Orbiter": ("Solar Orbiter", "2020-02-11", "tab:orange", "D"),
    "STEREO-A": ("STEREO-A", "2006-10-27", "tab:green", "s"),
    "SDO": ("SDO", "2010-02-12", "black", "^"),
    "SOHO": ("SOHO", "1995-12-03", "tab:cyan", "v"),
    "BepiColombo": ("BepiColombo", "2018-10-21", "tab:brown", "P"),
}
SHOWN = ["Sun", "Earth", "Solar Orbiter", "PSP"]   # ticked at the start
BLOCK_DAYS = 120        # trajectories are downloaded in blocks of this many days
SAMPLE = "2h"           # trajectory sampling, interpolated in between
OUT_DIR = Path(__file__).resolve().parent
POSITIONS_JSON = OUT_DIR / "cache" / "horizons_positions.json"
REFRESH_PREDICTIONS = 30 * u.day
HEADER = f"{'':<14}{'r (AU)':>7}{'r (Rsun)':>10}{'lon (deg)':>11}{'lat':>8}"
FOOTER = "lon = angle from Earth (+ west)"
RMAX = 1.75             # AU, radial extent of the plot (Mars reaches 1.67)
R_SS = 2.5 * u.R_sun    # PFSS source surface: the Parker spirals start here
V_RANGE = (300, 800)    # km/s, default solar wind speed range of the Parker spirals / connectivity


def spiral_shift(r_out, r_in, v):
    """Longitude (deg) by which a Parker spiral field line, for wind speed v (km/s), lies further
    west at r_in than at r_out (r_out > r_in): the Sun turns under the radially moving plasma."""
    return (sidereal_rotation_rate * (r_out - r_in) / (v * u.km / u.s)).to_value(u.deg)


def parse_speeds(text):
    """'300-800', '300 800' or '450' (km/s) -> (slow, fast); ValueError if not readable."""
    values = sorted(float(v) for v in text.replace("-", " ").replace(",", " ").split())
    if not 1 <= len(values) <= 2 or values[0] <= 0:
        raise ValueError(text)
    return values[0], values[-1]


class Trajectories:
    """Spacecraft trajectories from JPL Horizons, kept as heliocentric inertial (non-rotating)
    x, y, z so that they can be interpolated linearly in time."""

    def __init__(self):
        self.tracks = {name: [] for name in SPACECRAFT}   # name: [(jd array, xyz array 3xN)]
        self.tried = set()                                 # (name, block) already loaded or failed
        self.missing = {name: [] for name in SPACECRAFT}  # name: [(jd0, jd1)] Horizons had no data for
        self.store = self._read_store()                    # what horizons_positions.json holds

    @staticmethod
    def _read_store():
        """{spacecraft: {block: {"fetched", "jd", "x_au", "y_au", "z_au"}}} from POSITIONS_JSON."""
        try:
            return json.loads(POSITIONS_JSON.read_text())["spacecraft"]
        except FileNotFoundError:
            return {}
        except (ValueError, KeyError) as e:
            print(f"Ignoring unreadable {POSITIONS_JSON.name}: {e}")
            return {}

    def _stored(self, name, k):
        block = self.store.get(name, {}).get(str(k))
        if block is None:
            return None
        fetched = Time(block["fetched"])
        stop = Time((k + 1) * BLOCK_DAYS, format="jd")
        if stop > fetched and Time.now() - fetched > REFRESH_PREDICTIONS:   # an old prediction
            return None
        return np.array(block["jd"]), np.array([block["x_au"], block["y_au"], block["z_au"]])

    def _save(self, name, k, track):
        store = self._read_store()   # another window may have added blocks meanwhile
        jd, xyz = track
        store.setdefault(name, {})[str(k)] = {
            "fetched": Time.now().isot[:19],
            "jd": np.round(jd, 8).tolist(),
            **{f"{c}_au": np.round(v, 10).tolist() for c, v in zip("xyz", xyz)},
        }
        content = {"about": "JPL Horizons positions saved by plot_positions.py: per spacecraft and "
                            f"{BLOCK_DAYS}-day block (JD // {BLOCK_DAYS}), heliocentric inertial x, y, z (AU) "
                            f"every {SAMPLE} at the times jd (UTC Julian dates)",
                   "spacecraft": store}
        tmp = POSITIONS_JSON.with_suffix(".tmp")
        tmp.write_text(json.dumps(content))
        tmp.replace(POSITIONS_JSON)
        self.store = store

    @staticmethod
    def _query(horizons_name, start, stop):
        for attempt in range(2):   # Horizons sometimes answers an error page; try once more
            try:
                c = get_horizons_coord(horizons_name, {"start": start, "stop": stop, "step": SAMPLE})
            except Exception:
                if attempt == 0:
                    time.sleep(2)
                continue
            hci = c.transform_to(HeliocentricInertial(obstime=c.obstime))
            return c.obstime.jd, hci.cartesian.xyz.to_value(u.AU)
        return None

    def _load_block(self, name, k):
        """Trajectory for block k (days [k, k+1) * BLOCK_DAYS in JD), from the json or Horizons."""
        track = self._stored(name, k)
        if track is not None:
            return track
        horizons_name, first, _, _ = SPACECRAFT[name]
        start = max(Time(k * BLOCK_DAYS, format="jd"), Time(first))
        stop = Time((k + 1) * BLOCK_DAYS, format="jd") + 2 * u.h   # overlap the next block
        track = self._query(horizons_name, start, stop)
        if track is not None:
            self._save(name, k, track)
        return track

    def _covered(self, name, jd):
        for track in self.tracks[name]:
            if track[0][0] <= jd <= track[0][-1]:
                return track
        return None

    def fetch(self, t, names, status=None):
        """Load what spacecraft names need for time t (one spacecraft after the other)."""
        k = math.floor(t.jd / BLOCK_DAYS)
        for name in names:
            horizons_name, first, _, _ = SPACECRAFT[name]
            if (t < Time(first) or self._covered(name, t.jd) is not None
                    or any(j0 <= t.jd <= j1 for j0, j1 in self.missing[name])):
                continue
            if status:
                saved = (name, k) not in self.tried and self._stored(name, k) is not None
                status(f"Loading {name} from {POSITIONS_JSON.name if saved else 'JPL Horizons'} ...")
            track = None
            if (name, k) not in self.tried:
                self.tried.add((name, k))
                track = self._load_block(name, k)
            if track is None:   # e.g. near the end of the predictions: try just around t
                track = self._query(horizons_name, max(t - 3 * u.day, Time(first)), t + 3 * u.day)
            if track is None or not track[0][0] <= t.jd <= track[0][-1]:
                self.missing[name].append(((t - 1 * u.day).jd, (t + 1 * u.day).jd))
                print(f"Horizons: no {name} position for {t.iso[:16]}")
            else:
                self.tracks[name].append(track)

    def position(self, name, t):
        """Stonyhurst coordinate of spacecraft name at t, or None."""
        track = self._covered(name, t.jd)
        if track is None:
            return None
        jd, xyz = track
        x, y, z = (np.interp(t.jd, jd, xyz[k]) for k in range(3))
        hci = HeliocentricInertial(CartesianRepresentation(x, y, z, unit=u.AU), obstime=t)
        return hci.transform_to(HeliographicStonyhurst(obstime=t))


class PositionsPlot:
    def __init__(self, t, step):
        self.t = t
        self.traj = Trajectories()
        self.style = ({"Sun": ("orange", "o")} | {n: v[1:] for n, v in PLANETS.items()}
                      | {n: v[2:] for n, v in SPACECRAFT.items()})   # name: (colour, marker)
        self.names = list(self.style)

        self.fig = plt.figure(figsize=(12.5, 8.5))
        self.fig.canvas.manager.set_window_title("Planet and spacecraft positions")
        self.ax = self.fig.add_axes([0.03, 0.2, 0.52, 0.72], projection="polar")
        self.status = self.fig.text(0.03, 0.14, "", fontsize=9, color="0.3")
        self.v_range = V_RANGE

        # the table: one row per object, with a tick box to show it on the plot
        tax = self.fig.add_axes([0.51, 0.5, 0.45, 0.44], frameon=False)
        colors = [c for c, _ in self.style.values()]
        self.table = CheckButtons(tax, self.names, [n in SHOWN for n in self.names],
                                  label_props={"family": ["monospace"] * len(colors), "fontsize": [9.5] * len(colors),
                                               "color": colors},
                                  check_props={"facecolor": colors, "linewidth": 2})
        self.table.on_clicked(lambda _: self.update(fetch=True))
        mono = dict(transform=tax.transAxes, family="monospace", fontsize=9.5, va="center")
        tax.text(0.25, 1.0, HEADER, **mono)
        tax.text(0.25, 0.0, FOOTER, **mono, color="0.4")

        self.date_box = TextBox(self.fig.add_axes([0.12, 0.08, 0.22, 0.045]), "date (UTC) ")
        self.date_box.on_submit(self._date_typed)
        self.buttons = []
        for label, x0, action in [("<<", 0.36, lambda _: self.step(-1)), (">>", 0.42, lambda _: self.step(+1)),
                                  ("now", 0.48, lambda _: self.go(Time.now())),
                                  ("save png", 0.57, lambda _: self.save("png")),
                                  ("save pdf", 0.655, lambda _: self.save("pdf"))]:
            b = Button(self.fig.add_axes([x0, 0.08, 0.05 if len(label) < 4 else 0.075, 0.045]), label)
            b.on_clicked(action)
            self.buttons.append(b)
        self.spiral_toggle = CheckButtons(self.fig.add_axes([0.03, 0.005, 0.16, 0.05], frameon=False),
                                          ["Parker spirals"], [True])
        self.spiral_toggle.on_clicked(lambda _: self.update())
        self.speed_box = TextBox(self.fig.add_axes([0.3, 0.015, 0.1, 0.04]), "wind speed (km/s) ",
                                 initial=f"{V_RANGE[0]}-{V_RANGE[1]}")
        self.speed_box.on_submit(self._speeds_typed)
        self.fig.text(0.86, 0.32, "step", ha="left", va="bottom")
        self.step_buttons = RadioButtons(self.fig.add_axes([0.86, 0.02, 0.12, 0.3], frameon=False),
                                         list(STEPS), active=list(STEPS).index(step))
        self.go(t)

    def _speeds_typed(self, text):
        try:
            self.v_range = parse_speeds(text)
        except ValueError:
            self.set_status(f"Could not read the speeds {text!r} (e.g. 300-800)")
            return
        self.update()

    def set_status(self, text):
        self.status.set_text(text)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def shown(self):
        return [n for n, on in zip(self.names, self.table.get_status()) if on]

    def step(self, sign):
        self.go(self.t + sign * STEPS[self.step_buttons.value_selected])

    def _date_typed(self, text):
        try:
            t = parse_time(text.strip())
        except Exception:
            self.set_status(f"Could not read the date {text!r} (e.g. 2022-12-11 13:16)")
            return
        if t.isot[:16] != self.t.isot[:16]:
            self.go(t)

    def go(self, t):
        self.t = Time(t, scale="utc")
        self.update(fetch=True)

    def update(self, fetch=False):
        shown = self.shown()
        if fetch:
            self.traj.fetch(self.t, [n for n in shown if n in SPACECRAFT], status=self.set_status)
        coords = {name: get_body_heliographic_stonyhurst(body, time=self.t)
                  for name, (body, _, _) in PLANETS.items()}
        coords |= {name: self.traj.position(name, self.t) for name in SPACECRAFT}
        self.coords = coords
        self.plot(coords, shown)

    def draw(self, ax, coords, shown, rmax=None):
        """The polar plot of the shown objects; rmax None: the fixed RMAX of the window."""
        ax.clear()
        if self.spiral_toggle.get_status()[0]:
            for name, c in coords.items():
                if c is not None and name in shown and (name in SPACECRAFT or name == "Earth"):
                    self.draw_spirals(ax, c, self.style[name][0])
        if "Sun" in shown:
            ax.plot(0, 0, "o", color="orange", ms=12, label="Sun")
        for name, c in coords.items():
            if c is None or name not in shown:
                continue
            lon, r = c.lon.to_value(u.rad), c.radius.to_value(u.AU)
            color, marker = self.style[name]
            if name in SPACECRAFT:
                ax.plot([0, lon], [0, r], color=color, lw=0.8, alpha=0.6)
            ax.plot(lon, r, marker, color=color, ms=10 if name in PLANETS else 6, label=name,
                    mec="white" if name in SPACECRAFT else color, mew=0.5)
        if rmax is None:
            ax.set_rmax(RMAX)
            ax.set_rticks([0.25, 0.5, 0.75, 1, 1.25, 1.5])
        else:
            ax.set_rmax(rmax)
        ax.set_rlabel_position(-22.5)
        ax.set_thetagrids(range(0, 360, 45), [f"{a if a <= 180 else a - 360}°" for a in range(0, 360, 45)])
        ax.set_title(f"Stonyhurst heliographic   {self.t.iso[:16]} UTC", pad=18)

    def draw_spirals(self, ax, c, color):
        """Parker spirals from c down to R_SS for the slowest and fastest speed, band shaded."""
        r_sc = c.radius
        r = np.linspace(R_SS.to_value(u.AU), r_sc.to_value(u.AU), 200) * u.AU
        theta = [np.radians(c.lon.to_value(u.deg) + spiral_shift(r_sc, r, v)) for v in self.v_range]
        ax.fill(np.r_[theta[0], theta[1][::-1]], np.r_[r.value, r.value[::-1]], color=color, alpha=0.12, lw=0)
        for th, ls in zip(theta, ("--", ":")):
            ax.plot(th, r.value, ls, color=color, lw=0.9)

    @staticmethod
    def row(name, c):
        r = c.radius
        return (f"{name:<14}{r.to_value(u.AU):7.3f}{r.to_value(u.R_sun):10.1f}"
                f"{c.lon.wrap_at(180 * u.deg).to_value(u.deg):+11.2f}{c.lat.to_value(u.deg):+8.2f}")

    def plot(self, coords, shown):
        self.draw(self.ax, coords, shown)
        if shown:
            self.ax.legend(loc="lower left", bbox_to_anchor=(1.0, -0.02), fontsize=9, frameon=False, ncol=2)
        for name, label in zip(self.names, self.table.labels):
            c = coords.get(name)
            if name == "Sun":
                text = f"{name:<14}"
            elif c is not None:
                text = self.row(name, c)
            elif name in shown:
                text = f"{name:<14}{'not available':>20}"
            else:
                text = f"{name:<14}{'(tick to load)':>21}"
            label.set_text(text)
        self.date_box.eventson = False
        self.date_box.set_val(self.t.iso[:16])
        self.date_box.eventson = True
        self.set_status("")

    def save(self, fmt):
        """A figure of the ticked objects only, scaled to the farthest of them."""
        shown = [n for n in self.shown() if n == "Sun" or self.coords.get(n) is not None]
        rows = [n for n in shown if n != "Sun"]
        far = max((self.coords[n].radius.to_value(u.AU) for n in rows), default=0)
        fig = Figure(figsize=(11, 6))
        ax = fig.add_axes([0.03, 0.06, 0.5, 0.82], projection="polar")
        self.draw(ax, self.coords, shown, rmax=1.1 * far if far > 0 else RMAX)
        ax.text(np.deg2rad(-22.5), ax.get_rmax() * 1.13, "AU", ha="center", va="center", fontsize=9)
        y = 0.86
        fig.text(0.56, y, HEADER, family="monospace", fontsize=10)
        for name in rows:
            y -= 0.045
            fig.text(0.56, y, self.row(name, self.coords[name]), family="monospace", fontsize=10,
                     color=self.style[name][0])
        fig.text(0.56, y - 0.06, FOOTER, family="monospace", fontsize=9, color="0.4")
        if self.spiral_toggle.get_status()[0]:
            y -= 0.04
            fig.text(0.56, y - 0.06, f"shaded: Parker spirals for {self.v_range[0]:g}-{self.v_range[1]:g} km/s "
                                     f"(dashed / dotted) down to 2.5 Rsun", family="monospace", fontsize=9, color="0.4")
        if shown:
            ax.legend(loc="upper left", bbox_to_anchor=(1.05, (y - 0.1 - 0.06) / 0.82), fontsize=10,
                      frameon=False, ncol=2)
        out = OUT_DIR / f"positions_{self.t.strftime('%Y%m%d_%H%M')}.{fmt}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"Saved {out}")
        self.set_status(f"Saved {out.name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--date", help="start date (UTC), e.g. '2022-12-11 13:16' (default: now)")
    parser.add_argument("--step", default="1d", choices=list(STEPS), help="initial step (default 1d)")
    args = parser.parse_args()
    t = parse_time(args.date) if args.date else Time.now()
    plot = PositionsPlot(t, args.step)  # keep a reference: matplotlib holds callbacks only weakly
    plt.show()


if __name__ == "__main__":
    main()
