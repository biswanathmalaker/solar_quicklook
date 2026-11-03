#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
Download AIA full-disk files (level 1, aia.lev1_euv_12s) from JSOC via Fido.
Only the image segment is requested (JSOC would otherwise also send a "spikes"
file per frame).

Files go to {name}/{wavelength}/ next to this script, named
{name}_{wavelength}_YYYYMMDD_HHMMSS.fits (DATE-OBS); plot_ref_file.py shows
that folder and saves its figures there too.

  --ref-date only    one frame: the one closest to ref_date in [ref_date, ref_date + 12 s]
  --end-date         one frame every --cadence (default 1h) from ref_date to end_date,
                     in a single JSOC export per wavelength
  --wavelength       one or more of 171 (default), 193, 211, 304, 335, 94, 131

Frames already in the folder (within 60 s, or half the cadence, of a wanted time)
are not downloaded again unless --overwrite. A full-disk file is about 11 MB.

Uses the light_pipeline venv (shebang) and its mails.py.

Run:  ./download_ref_file.py --name plume_002 --ref-date 2019-05-02T18:00:00
      ./download_ref_file.py --name eq-ch --ref-date 2025-03-23T06:00 --end-date 2025-03-24T12:00 --wavelength 171 193
      ./download_ref_file.py --name eq-ch --ref-date 2025-03-23T06:00 --end-date 2025-03-24T12:00 --cadence 30min
"""
import argparse
import math
import re
import shutil
import socket
import sys
import tempfile
from pathlib import Path

import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from sunpy.net import Fido
from sunpy.net import attrs as a

LIGHT_PIPELINE = Path.home() / "projects/project3/light_pipeline"
sys.path.insert(0, str(LIGHT_PIPELINE))
from mails import mails  # noqa: E402  (same list as Functions.mails)

OUT_DIR = Path(__file__).resolve().parent
EMAIL = mails[3]
WINDOW = 12 * u.s
WAVELENGTHS = [171, 193, 211, 304, 335, 94, 131]
TIME_IN_NAME = re.compile(r"_(\d{8}_\d{6})\.fits$")


def data_dir(name, wavelength):
    return OUT_DIR / name / str(wavelength)


def obs_time(path):
    with fits.open(path) as hdul:
        header = hdul[-1].header  # JSOC files are compressed: the image header is in the last HDU
    return Time(header["DATE-OBS"])


def existing_times(folder):
    """DATE-OBS of the frames already in folder (from the file names)."""
    times = []
    for f in folder.glob("*.fits"):
        match = TIME_IN_NAME.search(f.name)
        times.append(Time.strptime(match.group(1), "%Y%m%d_%H%M%S") if match else obs_time(f))
    return times


def search(t0, t1, wavelength, cadence=None):
    attrs = [a.Time(t0, t1), a.jsoc.Series("aia.lev1_euv_12s"), a.Wavelength(wavelength * u.AA),
             a.jsoc.Segment("image"), a.jsoc.Notify(EMAIL)]
    if cadence is not None:
        attrs.append(a.Sample(cadence))
    result = Fido.search(*attrs)
    return result, (len(result[0]) if len(result) else 0)


def download(name, wavelength, ref_date, end_date=None, cadence=1 * u.h, overwrite=False, folder=None):
    """folder: where to put the frames (default {name}/{wavelength}/ here); returns the
    paths of the frames saved."""
    folder = Path(folder) if folder else data_dir(name, wavelength)
    folder.mkdir(parents=True, exist_ok=True)
    t0 = Time(ref_date)
    if end_date is None:
        wanted, tol = [t0], WINDOW
    else:
        t1 = Time(end_date)
        if t1 <= t0:
            raise SystemExit("--end-date must be after --ref-date")
        n = math.floor(((t1 - t0) / cadence).decompose().value + 1e-9) + 1
        wanted, tol = [t0 + k * cadence for k in range(n)], min(cadence / 2, 60 * u.s)
    have = [] if overwrite else existing_times(folder)
    missing = [t for t in wanted if not any(abs((h - t).to(u.s)) <= tol for h in have)]
    if not missing:
        print(f"{wavelength} Å: all {len(wanted)} frame(s) already in {folder}")
        return []
    print(f"{wavelength} Å: {len(missing)} of {len(wanted)} frame(s) to download "
          f"(about {11 * len(missing)} MB) -> {folder}")

    if end_date is None:
        result, n = search(t0, t0 + WINDOW, wavelength)
    else:   # one export over the missing span (frames in it that we have are skipped below)
        result, n = search(missing[0], missing[-1] + WINDOW, wavelength, cadence)
    if n == 0:
        print(f"{wavelength} Å: no frame found at JSOC for {missing[0].isot} - {missing[-1].isot}")
        return []

    saved = []
    with tempfile.TemporaryDirectory(dir=OUT_DIR, prefix=f".{name}_{wavelength}_") as tmp:
        files = Fido.fetch(result, path=str(Path(tmp) / "{file}"))
        if files.errors:
            print(f"{wavelength} Å: some downloads failed: {files.errors}")
        if end_date is None and files:
            files = [min(files, key=lambda f: abs((obs_time(f) - t0).sec))]
        for f in files:
            t = obs_time(f)
            if any(abs((h - t).to(u.s)) <= tol for h in have):
                continue
            path = folder / f"{name}_{wavelength}_{t.strftime('%Y%m%d_%H%M%S')}.fits"
            shutil.move(f, path)
            saved.append(path)
    print(f"{wavelength} Å: saved {len(saved)} frame(s) in {folder}")
    return saved


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--name", required=True, help="e.g. plume_002: the folder {name}/{wavelength}/")
    parser.add_argument("--ref-date", required=True, help="e.g. 2019-05-02T18:00:00 (UTC)")
    parser.add_argument("--end-date", help="download a frame every --cadence from --ref-date to this date (UTC)")
    parser.add_argument("--cadence", default="1h", help="with --end-date: e.g. 1h (default), 30min, 12min, 2h")
    parser.add_argument("--wavelength", type=int, nargs="+", default=[171], choices=WAVELENGTHS,
                        help="AIA channel(s) in Å (default 171), e.g. --wavelength 171 193")
    parser.add_argument("--overwrite", action="store_true", help="download again frames already there")
    args = parser.parse_args()
    try:
        cadence = u.Quantity(args.cadence).to(u.s)
    except (TypeError, ValueError, u.UnitsError):
        raise SystemExit(f"--cadence {args.cadence!r}: use e.g. 1h, 30min, 12min (minutes are 'min')")
    if cadence < 12 * u.s:
        raise SystemExit("--cadence must be at least 12 s (the AIA EUV cadence)")
    socket.setdefaulttimeout(300)  # JSOC can be slow to answer
    for wavelength in dict.fromkeys(args.wavelength):   # in the given order, without repeats
        download(args.name, wavelength, args.ref_date, args.end_date, cadence, args.overwrite)


if __name__ == "__main__":
    main()
