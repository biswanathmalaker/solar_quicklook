# Commands

Run from `~/projects/solar_quicklook`. Times are UTC.

## AIA frames and plume coordinates

| Command | What it does |
|---|---|
| `./download_ref_file.py --name plume_002 --ref-date 2019-05-02T18:00:00` | Download the AIA 171 frame at that time → `plume_002/171/` |
| `... --wavelength 171 193` | Download these channels (171 193 211 304 335 94 131) |
| `... --end-date 2019-05-03T18:00:00` | Download one frame per `--cadence` up to this date |
| `... --cadence 30min` | Cadence with `--end-date` (default 1h) |
| `... --overwrite` | Download again frames already there |
| `./plot_ref_file.py --name plume_002` | Show `plume_002/171/`; click the plume, and its coordinate goes to `ref_coords.csv` |
| `./plot_ref_file.py --name plume_002 --wavelength 193` | Same, for `plume_002/193/` |

Keys in `plot_ref_file.py`: left click = mark, `u` = undo, `b` / `n` = previous / next frame.

## Spacecraft positions

| Command | What it does |
|---|---|
| `./plot_positions.py` | Polar plot of planets and spacecraft now |
| `./plot_positions.py --date "2025-10-20 12:00" --step 1d` | Same at a date, stepping by 1h, 4h, 12h, 1d, 3d, 1wk or 1mo |

## Browse through time

| Command | What it does |
|---|---|
| `./browse_web.py` | Web browser at http://localhost:5006/browse, last date shown (settings kept) |
| `./browse_web.py 2025-10-20T12:00` | Same, starting at this date |
| `... --port 5007 --no-show` | Other port; do not open a browser tab |
| `./browse.py` | Old matplotlib browser (being replaced) |

In the web page: timeline (hover = details, click = go there, wheel = zoom time, line below = cursor time and its distance from the shown time); *display of* = per image (AIA channel, FSI, each original): stretch (linear, sqrt, square, log, asinh, sinh, power, histogram eq.), vmin / vmax (blank = auto: the central % of values), colour map, reversed, *default*; *search metadata* box (third row) + instrument buttons + Enter = keywords in all cached metadata, comma = OR, two or more instruments = only times all observed together (AND); results in the *search* tab, click a result to go to its start; *common FOV of* (fourth row) + instruments, *without* + instruments that must not observe that FOV (red, struck through; their times are cut out; one instrument is then enough), from / to (blank = all), min % (common area), min μ, *IRIS mosaics* + *find common FOV* = sets of observations at the same time on a common part of the Sun, results in the *common FOV* tab, click to go there; the four control rows stay at the top when scrolling; images: wheel = zoom, drag = pan, hover inside a FOV = its block, line below = cursor helioprojective x, y and Stonyhurst lon, lat, μ; date box + Enter; ◀ / ▶ = step by the cadence; IRIS slit-jaw (pale) and raster (coral) boxes; tabs: summary, Solar Orbiter files, DKIST, IRIS (programme, channels, level-2 links), search, positions; *positions* = PSP / Solar Orbiter sub-points, limbs and source regions + polar plot; *download original data* = full-resolution AIA (JSOC, else the VSO NASA mirror; if both fail the status line says so and the image stays synoptic) + Solar Orbiter files of the step; *show Solar Orbiter originals* = HRIEUV / HRI-Lya / PHI-HRT files nearest the time, full resolution.

In the old window: date box + Enter = go to a date; `<<` / `>>`, `b` / `n` or ← / → = step by the cadence; cadence and AIA channel buttons; *show* = which FOVs to draw; *positions* = spacecraft positions window; *original data* = download full-resolution AIA and the Solar Orbiter files of the step.

## Solar Orbiter: what was observed

| Command | What it does |
|---|---|
| `./find_obs.py 2025-10-01 2025-10-31` | HRIEUV and PHI-HRT observing blocks: time, cadence, SOOP, pointing, latitudes covered, Earth visibility |
| `... --products hri174 hriopn hrilya phi-hrt fsi174 fsi304 phi-fdt spice` | Choose the products (`hriopn` = HRIEUV 174 Å, filter wheel open) |
| `... --soop Polar` | Only blocks whose SOOP name contains "Polar" |
| `... --min-lat 60` | Only pointed blocks reaching \|lat\| ≥ 60° |
| `... --min-files 10` | Only blocks with at least 10 files |
| `... --both hri174 phi-hrt` | Only blocks where the two overlap in time and on the Sun |
| `... --csv blocks.csv` | Also save the table |
| `... --refresh` | Query SOAR again (new releases) instead of the cache |
| `... --no-headers` | Skip the slow FITS-header queries |
| `... --no-soop-info` | Do not list the SOOP summary entries |
| `./soar_meta.py --all` | Fetch the metadata the browsers use for every month since 2020-02 not yet in the cache (hours the first time; run again for months that failed) |
| `./soar_meta.py --from 2026-05 --to 2026-11 --refresh` | Query these months again (newly released data) |

## Fields of view on the full disk

| Command | What it does |
|---|---|
| `./plot_fov.py 2025-10-20T12:30` | HRIEUV / PHI-HRT / DKIST FOVs, PSP and its source region, on AIA 193 and FSI 174 → `fov/fov_20251020_1230.png` |
| `... --products hri174 hrilya phi-hrt` | Which Solar Orbiter FOVs to draw |
| `... --window 24h` | Use Solar Orbiter files within this of the time (default 2h) |
| `... --aia 171` | AIA channel (default 193) |
| `... --fsi 304` | FSI channel (default 174) |
| `... --no-fsi` | AIA panel only |
| `... --spacecraft PSP STEREO-A BepiColombo` | Mark these: sub-spacecraft point and limb |
| `... --connect PSP "Solar Orbiter" "Earth / L1"` | Draw their solar source regions (Parker spiral) |
| `... --speeds 350-750` | Wind speeds for the source regions (default 300-800 km/s) |
| `... --margin 5` | Degrees around the source regions (default 10) |
| `... --dkist-window 6h` | DKIST datasets within this of the time (default 3h) |
| `... --no-dkist` | Do not draw DKIST |
| `... --positions` | Add the spacecraft-positions panel |
| `... --zoom` | Zoom on the Solar Orbiter FOVs |
| `... --no-show` | Only save the figure |

## DKIST, IRIS, SOOPs, SPICE studies and keyword search

| Command | What it does |
|---|---|
| `./spice_meta.py solo_L2_spice-n-ras_20251020T113040_V05_352321876-008-DR6.fits` | The file's SPICE study: description, raster, exposure, spectral lines |
| `./dkist_meta.py 2025-10-01 2025-10-31` | DKIST datasets: instrument, times, FOV box, proposal, lines |
| `./dkist_meta.py 2025-10-01 2025-10-31 --refresh` | Same, querying DKIST again instead of the cache |
| `./dkist_meta.py --all` | Fetch every month since 2020-01 not yet in the cache (~2 min) |
| `./iris_meta.py 2025-10-10 2025-10-17` | IRIS observations: times, OBSID, pointing, roll, raster size, aim, planner |
| `... --pole` | Only polar pointings (\|y\| > 800″) |
| `... --long` | Also the programme, raster (steps, cadence, exposure), slit-jaw channels, data size and link |
| `... --refresh` | Query LMSAL HCR again instead of the cache |
| `./iris_meta.py --all` | Fetch every month since 2013-07 not yet in the cache (~40 s a month; run again for months that failed) |
| `./search_meta.py "coronal hole, plume"` | Observations in all cached metadata matching any of the comma-separated keywords (OR; case-insensitive; several words = that phrase) |
| `... --inst IRIS SPICE` | Only these (HRIEUV HRI-Lya PHI-HRT SPICE FSI IRIS DKIST); two or more = only the times when all of them have a matching observation at once (AND) |
| `... --from 2023-01 --to 2024-12` | Only this time range |
| `... --max 200` | Print this many (default 50; 0 = all) |
| `./fov_search.py PHI-HRT IRIS --from 2025-10 --to 2025-11` | Common field of view: sets of observations of all these instruments at the same time overlapping on the Sun (common time, common area as % of the smallest FOV, its centre) |
| `... --min 50` | Only common areas of at least 50 % of the smallest FOV |
| `... --mu 0.3` | Solar Orbiter vs Earth: common area seen at μ ≥ 0.3 from both (default 0.1; polar regions ~0.2-0.35 from Earth) |
| `... --mosaics` | Include IRIS full-disk mosaics (whole disk in HCR for hours, one tile at a time; left out by default) |
| `./fov_search.py HRIEUV --without IRIS` | Instruments that must NOT observe the common FOV: the times they share any part of it are cut out, pieces of at least 5 min kept; one instrument is then enough (here: HRIEUV with no IRIS on its FOV) |
| `./fov_search.py HRIEUV PHI-HRT --without IRIS DKIST --left 10` | Several on each side; `--left` = shortest piece kept, in minutes |
| `./fov_search.py --build` | (Re)build its index from the cache (`cache/fov_search/`, ~3 min the first time) |
| `./search_meta.py --build` | (Re)build the search index from the cache (`cache/search/`, ~2 min the first time) |
| `./soops.py` | All runs in ESA's SOOP summary |
| `./soops.py Polar` | Runs whose SOOP name contains "Polar" |
| `./soops.py --date 2025-10-20` | Runs covering this date |
| `./soops.py --pages` | Fetch the ESA description page of every SOOP (summary + SOAR names) not yet cached (`--refresh`: all again) |

## Update the cache with new releases

All metadata is in `cache/` (SOAR since 2020-02, DKIST since 2020-01, IRIS since 2013-07) and is never queried again by itself. Recent months were saved as they were (on 2026-11-04: SOAR 2026-05 to 2026-11 incomplete, not yet released; IRIS and DKIST 2026-11 partial). To update them, then restart `browse_web.py`:

| Command | What it does |
|---|---|
| `./soar_meta.py --from 2026-05 --to 2026-11 --refresh` | Solar Orbiter: all products the browsers use, FSI, SPICE studies |
| `./iris_meta.py 2026-10-01 2026-11-30 --refresh` | IRIS (also the month before, if the start is within 2 days of it) |
| `./dkist_meta.py 2026-11-01 2026-11-30 --refresh` | DKIST (also the month before, if the start is within 1 day of it) |
| `./soar_meta.py --all`, `./iris_meta.py --all`, `./dkist_meta.py --all` | Add the months since the last bulk download (cached months are skipped) |

## Help

| Command | What it does |
|---|---|
| `./help.py` | Print the main commands |
| `./<script>.py -h` | All options of a script |
