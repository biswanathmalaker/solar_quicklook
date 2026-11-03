#!/usr/bin/env python3
"""Print the few commands needed for one plume. Run: python3 help.py"""

print("""\
Everything for a name goes in <name>/<wavelength>/ (frames, saved figures).

1. Download one AIA frame (171 Å) at the plume's reference time (UTC) -> plume_002/171/:
     ./download_ref_file.py --name plume_002 --ref-date 2019-05-02T18:00:00

2. Other / several channels (171 193 211 304 335 94 131) and several frames, one per hour
   (or --cadence 30min ...); frames already there are not downloaded again:
     ./download_ref_file.py --name eq-ch --ref-date 2025-03-23T06:00 --end-date 2025-03-24T12:00 --wavelength 171 193

3. Show <name>/<wavelength>/ and click the plume; on closing, the last point goes to ref_coords.csv:
     ./plot_ref_file.py --name plume_002
     ./plot_ref_file.py --name plume_002 --wavelength 193

   left click = mark point   u = undo   sliders = contrast (applied on mouse release)
   mouse over image = toolbar shows x, y (arcsec) and Stonyhurst lon, lat (deg)
   stretch buttons = linear (default) / sqrt / log / asinh / squared
   inset = then click 2 corners (or type x, y + Enter) -> zoomed inset joined by lines;
           drag it to move, drag a corner to resize
   spacecraft tick boxes = sub-spacecraft point (x) + its limb (solid = seen by AIA) for PSP,
                           Solar Orbiter (default), STEREO-A, BepiColombo; Earth / L1
   connectivity region = part of the Sun each ticked one can be connected to (Parker
                         spiral to 2.5 Rsun for 300-800 km/s, ±10°; edit km/s and ± boxes)
                         style button: outline (default) / shaded / dim outside
   < prev / next > (keys b / n) = step through the frames (when there are several)
   save png / save pdf = figure without widgets, next to the frame (same name .png / .pdf)

4. Download again frames already there:  add --overwrite to command 1 or 2

5. See the saved coordinates:
     cat ref_coords.csv

6. Planet / spacecraft positions (Earth, SDO, SOHO, PSP, Solar Orbiter, STEREO-A, BepiColombo):
     ./plot_positions.py --date "2022-12-11 13:16"
   date box + Enter = go to a date   << / >> = step back / forward   step = 1h ... 1mo
   tick boxes in the table = show / hide an object (default: Sun, Earth, Solar Orbiter, PSP)
   Parker spirals + wind speed box = field-line band of each shown spacecraft / Earth (300-800 km/s)
   save png / save pdf = ticked objects only, scaled to the farthest -> positions_YYYYMMDD_HHMM.png

7. What Solar Orbiter observed (no FITS download; SOAR metadata cached in cache/soar/):
     ./find_obs.py 2025-10-01 2025-10-31                        # HRIEUV + PHI-HRT observing blocks
     ./find_obs.py 2025-01-01 2026-03-31 --products phi-hrt --soop Polar
     ./find_obs.py 2025-01-01 2025-12-31 --both hri174 phi-hrt  # both instruments, same time and place
   products: hri174 hrilya phi-hrt fsi174 fsi304 phi-fdt spice    --min-lat 60   --csv out.csv   --refresh

8. Fields of view at a time, on AIA (Earth view) and EUI/FSI (Solar Orbiter view) -> fov/fov_YYYYMMDD_HHMM.png:
   HRIEUV / PHI-HRT (+ DKIST) FOVs, PSP position + limb, PSP source region (Parker spiral, 300-800 km/s)
     ./plot_fov.py 2025-10-20T12:30
     ./plot_fov.py 2025-03-16T12:00 --aia 171 --fsi 304 --zoom
     ./plot_fov.py 2025-10-20T12:30 --positions --connect PSP "Solar Orbiter" "Earth / L1" --speeds 350-750
   --spacecraft PSP STEREO-A ...   --dkist-window 6h   --no-dkist   --window 24h

9. What DKIST observed (bounding boxes, instruments, lines):
     ./dkist_meta.py 2025-10-01 2025-10-31

10. ESA SOOP summary (from soops-summary.html; save a newer copy of the web page over it to update):
     ./soops.py Polar               ./soops.py --date 2025-10-20

11. Browse through time in the web browser (http://localhost:5006/browse; settings kept):
     ./browse_web.py                ./browse_web.py 2025-10-20T12:00
   timeline: hover = details, click = go there, wheel = zoom time; images: wheel zoom, drag pan,
   hover inside a FOV = its block; below each image = cursor x, y (arcsec), lon, lat, mu; ◀ / ▶ = step; tabs: summary, Solar Orbiter files, DKIST, positions
   positions = PSP / Solar Orbiter sub-points, limbs, source regions + polar plot
   download original data = full-res AIA (JSOC, 1-3 min) + Solar Orbiter files of the step (cache/)
   show Solar Orbiter originals = HRIEUV / HRI-Lya / PHI-HRT files nearest the time, full resolution
   Old matplotlib version (AIA + what Solar Orbiter / DKIST observed in each step):
     ./browse.py                    ./browse.py 2025-10-20T12:00
   date box + Enter, << / >> (b / n, arrows), cadence 12s .. 1wk, AIA channel, show (FOVs),
   positions (2nd window), original data (full-res AIA + Solar Orbiter files, kept in cache/)

Everything queried / downloaded is kept in cache/ and not fetched again (--refresh in find_obs.py
and dkist_meta.py to get newly released data). See README.md, "The cache".
""")
