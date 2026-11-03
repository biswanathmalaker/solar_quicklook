#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python
"""
What a SPICE study observes: its description, raster or sit-and-stare, slit, exposure,
and its spectral windows (lines), with the rough formation temperature of each line.

SOAR's header table (v_spi_sc_fits) gives the study name of every file (soar_meta.
spice_modes) but only the first spectral window. The windows come from the FITS headers of
one file of the study, read with HTTP range requests from the SPICE data release at IAS
(https://spice.osups.universite-paris-saclay.fr/spice-data/, a few kB per window, not the
~90 MB file). A study's windows are fixed, so this is done once per study and kept in
cache/spice/studies.json.

    import spice_meta as sp
    s = sp.study("SCI_POLAR_SS_SL04_59.6S_FF", "solo_L2_spice-n-ras_20251020T113040_V05_352321876-008-DR6.fits")
    sp.lines_text(s)      # 'Mg IX 70.6 (1 MK) | Ne VIII 77.0, 78.0 (0.6 MK) | ...'

Run:  ./spice_meta.py solo_L2_spice-n-ras_20251020T113040_V05_352321876-008-DR6.fits
"""
import argparse
import json
import re

import requests
from astropy.io.fits import Header

from soar_meta import OUT_DIR

SERVER = "https://spice.osups.universite-paris-saclay.fr/spice-data"
CACHE = OUT_DIR / "cache" / "spice" / "studies.json"
KEYS = ["STUDY", "STUDYDES", "OBS_DESC", "STUDYTYP", "DATAPROD", "PURPOSE", "SLIT_WID", "XPOSURE", "TELAPSE",
        "NRASTERS", "NWIN", "TARGET", "AUTHOR", "SOOPNAME"]
# approximate temperature of maximum formation (log T, ionisation equilibrium), by ion
LOG_T = {"H I": 4.0, "C II": 4.6, "O II": 4.5, "C III": 4.8, "N III": 4.9, "O III": 4.9, "S IV": 5.0,
         "N IV": 5.1, "O IV": 5.2, "S V": 5.2, "S VI": 5.3, "O V": 5.4, "O VI": 5.5, "Ne VI": 5.6,
         "Ne VII": 5.7, "Ne VIII": 5.8, "Mg VIII": 5.9, "Mg IX": 6.0, "Mg X": 6.1, "Si XII": 6.3,
         "Fe XVIII": 6.8, "Fe XX": 7.0}
_memo = {}


def _get(session, url, start, n):
    r = session.get(url, headers={"Range": f"bytes={start}-{start + n - 1}"}, timeout=60)
    if r.status_code == 416:
        return b""
    r.raise_for_status()
    if r.status_code != 206:
        raise RuntimeError("the server ignored the byte range")
    return r.content


def fits_headers(url, stop=None, max_hdus=60):
    """Headers of the HDUs of a remote FITS file, read with byte ranges (header blocks
    only, the data are skipped). stop(header) -> True ends the reading after that HDU."""
    s, off, out = requests.Session(), 0, []
    for _ in range(max_hdus):
        buf = b""
        while True:
            chunk = _get(s, url, off + len(buf), 2880 * 4)
            if not chunk:
                return out
            k0 = len(buf) - len(buf) % 80
            buf += chunk
            end = next((k for k in range(k0, len(buf), 80) if buf[k:k + 8] == b"END     "), None)
            if end is not None:
                break
        h = Header.fromstring(buf[:end + 80].decode("ascii", errors="replace"))
        out.append(h)
        if stop is not None and stop(h):
            return out
        size = 0
        if h.get("NAXIS", 0):
            size = abs(h.get("BITPIX", 8)) // 8
            for a in range(1, h["NAXIS"] + 1):
                size *= h[f"NAXIS{a}"]
        size = h.get("GCOUNT", 1) * (h.get("PCOUNT", 0) + size)
        off += (end // 2880 + 1) * 2880 + -(-size // 2880) * 2880
    return out


def file_url(filename):
    """URL of a SPICE L2 file in the latest data release that has it (None if none)."""
    rel = requests.get(f"{SERVER}/metadata/latest-release.txt", timeout=30).text.strip()
    m = re.search(r"_(\d{4})(\d{2})(\d{2})T", filename)
    for r in dict.fromkeys([rel, "6.0", "5.0"]):
        url = f"{SERVER}/release-{r}/level2/{m.group(1)}/{m.group(2)}/{m.group(3)}/{filename}"
        if requests.head(url, timeout=30).status_code == 200:
            return url
        # the release may hold another version of the file: same name up to _V
        day = requests.get(url.rsplit("/", 1)[0] + "/", timeout=60)
        stem = re.sub(r"_V\d+.*$", "", filename)
        other = re.findall(rf'href="({re.escape(stem)}_V[^"]+\.fits)"', day.text) if day.ok else []
        if other:
            return url.rsplit("/", 1)[0] + "/" + sorted(other)[-1]
    return None


def _load():
    return json.loads(CACHE.read_text()) if CACHE.is_file() else {}


def study(name, filename, refresh=False):
    """Description and windows of a study (OBSERVATION_MODE / STUDY), from the cache or
    from the headers of filename (a file of that study; name=None: the file's study). dict with the KEYS (lower case),
    nx, ny (slit positions, pixels along the slit), cdelt1, cdelt2 (arcsec), and windows:
    [dict(name, wavemin, wavemax) in nm]. None if it cannot be read (tried again next
    session)."""
    if name is not None and name in _memo and not refresh:
        return _memo[name]
    cache = _load()
    if name is not None and name in cache and not refresh:
        _memo[name] = cache[name]
        return cache[name]
    try:
        url = file_url(filename)
        if url is None:
            raise RuntimeError("not in a SPICE data release")
        nwin = []

        def stop(h):   # the windows come first (the primary HDU is the first): stop after the last
            if "WAVEMIN" in h and h.get("EXTNAME"):
                nwin.append(h)
            return bool(nwin) and len(nwin) >= int(nwin[0].get("NWIN", 99))
        hs = fits_headers(url, stop)
    except Exception as e:
        print(f"  SPICE study {name or filename}: {e}")
        if name is not None:
            _memo[name] = None
        return None
    h0 = nwin[0] if nwin else hs[0]
    name = name or str(h0.get("STUDY"))
    out = {k.lower(): h0.get(k) for k in KEYS}
    out.update(file=url.rsplit("/", 1)[1], nx=h0.get("NAXIS1"), ny=h0.get("NAXIS2"),
               cdelt1=h0.get("CDELT1"), cdelt2=h0.get("CDELT2"),
               windows=[dict(name=str(h["EXTNAME"]), wavemin=float(h["WAVEMIN"]), wavemax=float(h["WAVEMAX"]))
                        for h in nwin])
    cache = _load()
    cache[name] = out
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    tmp = CACHE.with_suffix(".tmp")
    tmp.write_text(json.dumps(cache, indent=1))
    tmp.replace(CACHE)
    _memo[name] = out
    return out


def window_lines(name):
    """'S V 786 / O IV 787 (Merged)' -> [('S V', 78.6), ('O IV', 78.7)]; 'Ly Beta 1025' ->
    [('H I', 102.5)] (wavelengths in nm)."""
    out = []
    for part in name.split("/"):
        part = part.replace("Ly Beta", "H I").replace("Ly Gamma", "H I").replace("Ly beta", "H I")
        m = re.search(r"\b([A-Z][a-z]?)\s+([IVXL]+)\s+(\d{3,4}(?:\.\d+)?)", part)
        if m:
            out.append((f"{m.group(1)} {m.group(2)}", float(m.group(3)) / 10))
    return out


def lines_text(s):
    """The windows of a study as 'ion wavelength (T)', grouped by ion."""
    if not s or not s.get("windows"):
        return ""
    by_ion = {}
    for w in s["windows"]:
        lines = window_lines(w["name"]) or [(w["name"].split(" - ")[0], (w["wavemin"] + w["wavemax"]) / 2)]
        for ion, wl in lines:
            by_ion.setdefault(ion, []).append(wl)
    parts = []
    for ion, wls in by_ion.items():
        t = LOG_T.get(ion)
        temp = "" if t is None else f" ({10 ** t / 1e6:.2g} MK)"
        parts.append(f"{ion} {', '.join(f'{x:.1f}' for x in sorted(set(wls)))}{temp}")
    return " | ".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("filename", help="a SPICE L2 file name (its study is looked up)")
    ap.add_argument("--study", help="study name (default: the file's STUDY keyword)")
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()
    s = study(a.study, a.filename, a.refresh)
    if s:
        print(json.dumps(s, indent=1))
        print(lines_text(s))


if __name__ == "__main__":
    main()
