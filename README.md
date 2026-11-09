# solar_quicklook

Quick-look and data-search tools for solar observations: SDO/AIA frames and plume
coordinates (below), and what Solar Orbiter (EUI, PHI) observed and where its
high-resolution telescopes pointed, drawn on the full disk from both vantage points
(see *Solar Orbiter: search and fields of view*). Copied from `~/quick_aia` on
2026-11-03 (code only; the AIA frames stayed there).

Originally: movies of the selected plumes (`plume_002.mp4`, `plume_003.mp4`, ...) plus two
small scripts to get each plume's **reference coordinate**: download AIA
full-disk frames (171 Å by default) around the plume's reference time, then click on the plume to
read its position.

The start/end times and reference time come from watching the movie; the
scripts only need the reference time.

## Files

| File | What it is |
|---|---|
| `plume_XXX.mp4` | Movie of one selected plume. |
| `download_ref_file.py` | Downloads AIA level-1 full-disk frames (JSOC, via Fido) → `<name>/<wavelength>/`: one at `--ref-date`, or one per `--cadence` (default 1h) up to `--end-date`; one or more wavelengths. |
| `plot_ref_file.py` | Shows the frames of `<name>/<wavelength>/` (steps through them if several); click the plume to get its coordinate → `ref_coords.csv`. |
| `plot_positions.py` | Interactive polar plot (Stonyhurst) of the inner planets and of PSP, Solar Orbiter, STEREO-A, SDO, SOHO, BepiColombo at a date; step back / forward by 1h ... 1mo. |
| `positions_YYYYMMDD_HHMM.png/.pdf` | Figures saved by `plot_positions.py`. |
| `cache/horizons_positions.json` | Spacecraft positions downloaded from JPL Horizons (by `plot_positions.py` and `plot_ref_file.py`), read first next time; a period missing from it is downloaded and added. |
| `<name>/<wavelength>/` | e.g. `plume_002/171/`, `eq-ch/193/`: the AIA frames `<name>_<wavelength>_YYYYMMDD_HHMMSS.fits` (DATE-OBS), and the figures saved from them by `plot_ref_file.py` (same name, `.png` / `.pdf`). |
| `ref_coords.csv` | One row per plume: `plume_name, date_obs, xcen, ycen` (made by `plot_ref_file.py`). |
| `browse_web.py` | Web version of `browse.py` (Panel + Bokeh, http://localhost:5006/browse; settings in `browse_web_settings.json`): timeline of every observing block (HRIEUV, HRI-Lya, PHI-HRT, SPICE), FSI run and DKIST dataset (hover for details, click to go there; below it the time under the cursor and its distance from the shown time); a *display of* row: for each image on the page (AIA channel, FSI, each Solar Orbiter original) its stretch (linear, sqrt, square, log, asinh, sinh, power, histogram equalisation), vmin / vmax (blank = the central % of the values), colour map and reversed, kept per image kind in the settings (*default* = back to the reduced PNG's display); AIA and FSI with zoom / pan kept between steps and hover on each FOV for its block, and below each image the cursor position (helioprojective x, y; Stonyhurst lon, lat and μ for that image's observer, or the distance in R☉ off the disk); IRIS observations (timeline row; slit-jaw area pale, raster area coral, on both images); tabs: the summary of `browse.py` plus IRIS, the Solar Orbiter files and the DKIST datasets of the step, IRIS (aim, planner, programme, raster and slit-jaw channels, lines, level-2 links and sizes), search (results of the *search metadata* box in the third control row: keywords in all cached metadata of the chosen instruments, `search_meta.py`; click a result to go to its start), positions. A fourth row, *common FOV of* (instruments, *without* (instruments that must not observe that FOV), from / to, min % (common area), *find common FOV*), runs `fov_search.py`; results in the *common FOV* tab, click to go there. The four control rows (date and buttons, display, search, common FOV) stay at the top of the window when the page is scrolled. Dates before a mission's start (Solar Orbiter 2020-02, DKIST 2020-01, IRIS 2013-07) query nothing for it. Buttons: *positions* (PSP / Solar Orbiter sub-points, limbs and source regions on both images; polar plot with Parker spirals), *download original data* (full-resolution AIA from JSOC, else from the VSO NASA mirror (sdo7.nascom.nasa.gov); if both fail the status line says so and stays while that time is shown, then shown at 2048 px, and the Solar Orbiter files of the step), *show Solar Orbiter originals* (HRIEUV / HRI-Lya / PHI-HRT B_LOS files nearest the time, full resolution, own view). Images reduced once to 1024-px PNGs (`cache/web/`), FOV outlines in closed form, neighbouring steps prepared in the background: ~0.1 s a step once cached. |
| `browse.py` | (Old, matplotlib) Interactive: step through time (date box, cadence 12s ... 1wk, AIA channel; settings kept in `browse_settings.json`); AIA with the FOVs of what Solar Orbiter (HRIEUV, HRI-Lya, PHI-HRT, SPICE) and DKIST observed in each step (one colour per instrument), FSI when Solar Orbiter observed, a summary (what each instrument measures, campaign aim, time range, cadence, FOV; SPICE study and lines; DKIST instrument, lines, proposal); "positions" window (spacecraft, limbs, source regions); "original data" (full-resolution AIA + the Solar Orbiter files). |
| `soar_meta.py` | Solar Orbiter archive (SOAR) metadata: the L2 files of a time range, their pointing, and observing blocks. Used by the two scripts below; also importable. |
| `find_obs.py` | Table of Solar Orbiter observing blocks in a time range: where HRIEUV / PHI-HRT pointed, latitudes covered, whether Earth sees the same region; filters by SOOP, latitude, overlap of two instruments. |
| `plot_fov.py` | Figure at a time: HRIEUV / HRI-Lya / PHI-HRT and DKIST fields of view, PSP (and other spacecraft) positions and the solar source regions their wind can come from, on AIA (Earth view) and EUI/FSI (Solar Orbiter view); optional positions panel → `fov/fov_YYYYMMDD_HHMM.png`. |
| `spice_meta.py` | What a SPICE study observes: description, raster/slit/exposure, spectral windows (lines) with formation temperatures. Read once per study from the FITS headers of one file (HTTP byte ranges on the SPICE data release at IAS, a few kB), cached in `cache/spice/studies.json`. `./spice_meta.py <SPICE L2 file name>`. |
| `search_meta.py` | Keyword search (case-insensitive; keywords separated by commas, any of them (OR); several words in one keyword = that phrase; with two or more instruments only the times when all of them have a matching observation going on at once (AND), one result per common interval) in all the cached metadata: Solar Orbiter observing blocks (instrument, SOOP name, the SOOP summary's description, coordinators and notes, cached SOOP pages, SPICE study names), FSI runs (SOOP names), IRIS observations (OBSID, aim, planners, programme, target, description, lines) and DKIST datasets (proposal, experiment description, target types, lines). The index is built per month into `cache/search/` (~1.5 min the first time) and rebuilt for a month when its metadata is newer. Used by the *search* tab of `browse_web.py`; also a command: `./search_meta.py "coronal hole" --inst IRIS DKIST`. |
| `fov_search.py` | Common field of view, apart from the keyword search: sets of observations, one per instrument (two or more), going on at the same time whose fields of view overlap on the Sun; gives the common time, the common area as a share of the smallest FOV and its centre (Stonyhurst lat, lon). FOVs: Solar Orbiter blocks (block pointing), FSI 3.8° square, IRIS slit-jaw area, DKIST bounding box. Instruments on one platform (all Solar Orbiter, or IRIS + DKIST from Earth) are compared directly in that view, off-limb parts included; Solar Orbiter against Earth through their footprints on the surface, seen from above the footprints (off-limb parts count at the limb; coronal overlap between the two views is not counted), and only surface both see at μ ≥ 0.1 (`--mu`), so when Solar Orbiter is round the far side a strip both see edge-on at a limb is not a common FOV; each result gives the μ of its centre from each side; overlaps under 1 % of the smallest FOV (edges touching) are dropped; IRIS full-disk mosaics (HCR gives the whole disk for hours, IRIS sees one tile at a time) are left out unless `--mosaics`; area measured on a 64×64 grid. Index (time, observer, FOV corners) per source and month in `cache/fov_search/`. `--without INST ...`: instruments that must NOT observe the common FOV; the times when one shares any part of it are cut out of each result and the pieces left of at least 5 min (`--left`) are listed, so one instrument is enough (e.g. `./fov_search.py HRIEUV --without IRIS`: HRIEUV with no IRIS on its FOV). FOVs only: an IRIS box inside the HRIEUV FOV but away from a given feature still excludes that time. The *common FOV* row and tab of `browse_web.py`; also a command: `./fov_search.py PHI-HRT IRIS --from 2025-10 --to 2025-11 --min 50`. |
| `iris_meta.py` | IRIS observations (aim, planner, programme, pointing, roll, raster and slit-jaw FOVs, steps, cadences, exposures, channels, lines, level-2 links) in a time range, from the LMSAL Heliophysics Coverage Registry (the database behind iris.lmsal.com/search, no login); cached in `cache/iris/`. Also a command: `./iris_meta.py 2025-10-10 2025-10-17 --pole --long`. |
| `dkist_meta.py` | DKIST datasets (instrument, times, bounding box, proposal, experiment description, target types, lines) in a time range, from the DKIST Data Center API; cached in `cache/dkist/`. Also a command: `./dkist_meta.py 2025-10-01 2025-10-31`. |
| `soops.py` | ESA's SOOP summary (`soops-summary.html`, saved by hand from the web page) as a table, `cache/soops/soops_summary.csv`; `./soops.py Polar`, `./soops.py --date 2025-10-20`. `find_obs.py` lists the runs behind its blocks. `soops.page(soop)` reads a SOOP's ESA description page (aim, instrument modes), cached in `cache/soops/pages/`. |
| `soops-summary.html` | https://www.cosmos.esa.int/web/solar-orbiter/soops-summary saved from the browser; save a newer copy over it to update (re-parsed automatically). |
| `cache/` | Everything queried or downloaded, kept so it is never fetched twice (see *The cache* below). |

## Steps for one plume

```bash
cd ~/projects/solar_quicklook

# 1. watch plume_002.mp4 and note its reference time, then:
./download_ref_file.py --name plume_002 --ref-date 2019-05-02T18:00:00          # -> plume_002/171/
#    (other / several channels: --wavelength 171 193 -> plume_002/171/, plume_002/193/)
#    (several frames, one per hour: add --end-date 2019-05-03T18:00:00)

# 2. click on the plume (its base) to read the coordinate:
./plot_ref_file.py --name plume_002    # add --wavelength 193 for plume_002/193/
```

### 1. `download_ref_file.py`

- Searches JSOC (`aia.lev1_euv_12s`) for each `--wavelength` given (one or more
  of 171 (default), 193, 211, 304, 335, 94, 131 Å) and saves the frames in
  `<name>/<wavelength>/` as `<name>_<wavelength>_YYYYMMDD_HHMMSS.fits` (DATE-OBS).
- `--ref-date` only: one frame, the one closest to `--ref-date` in
  [`--ref-date`, `--ref-date` + 12 s].
- With `--end-date`: one frame every `--cadence` (default `1h`; e.g. `30min`,
  `12min`, `2h` -- minutes are `min`, `m` would be metres) from `--ref-date` to
  `--end-date`, in a single JSOC export per wavelength. A day at 1h is 24
  files, about 270 MB per wavelength.
- Frames already in the folder (within 60 s, or half the cadence, of a wanted
  time) are not downloaded again, so extending a range only fetches the new
  frames; `--overwrite` downloads them all again.
- Only the image segment is downloaded (JSOC otherwise also sends a "spikes"
  file per frame).
- Times are UTC. The JSOC notification email is `mails[3]`.
- Takes 1-3 minutes per wavelength (most of it waiting for JSOC to prepare the
  export). A full-disk file is about 11 MB.
- Run it as `./download_ref_file.py` (the light_pipeline Python), not with the
  system `python3` (see Environment below).

### 2. `plot_ref_file.py`

| Action | What it does |
|---|---|
| left click | Mark a point; its helioprojective x, y (arcsec, at the file's DATE-OBS) is printed and written on the image. The newest point is red, earlier ones cyan. |
| mouse over the image | The toolbar shows helioprojective `x`, `y` (arcsec) and heliographic Stonyhurst `lon`, `lat` (deg), or *off disk*; also over the inset. |
| `u` | Remove the last point. |
| vmin / vmax sliders | Image contrast (DN/s). The image is redrawn when you release the mouse (a full redraw takes ~1.5 s), so drag, then let go. |
| stretch buttons | Colour scale: linear (default), sqrt, log, asinh or squared. |
| inset | Then click two opposite corners of a region on the image: they fill the *bottom left* / *top right* boxes (arcsec) and a zoomed inset of the region appears (top right), with a box on the image and two lines joining the box to the inset. You can also type `x, y` in the boxes and press Enter. Drag the inset to move it; drag one of its corners to resize it (the shape is kept). The cursor shows which: a move cursor inside, a crosshair at the corners. *remove inset* removes it. The inset follows the sliders and stretch, and is in the saved figures. |
| spacecraft tick boxes (right) | PSP, Solar Orbiter (ticked at the start), STEREO-A, BepiColombo: draw the **sub-spacecraft point** (`x`, the centre of the disk that spacecraft sees; `o` at the limb when it is on the far side) and the **spacecraft's limb** (solid where AIA sees it, dotted where hidden): the region inside it is seen by both. The panel lists distance, angle from SDO and Stonyhurst lon / lat. Positions come from JPL Horizons at DATE-OBS (shared cache with `plot_positions.py`; ~4 s per spacecraft the first time). |
| connectivity region (right panel: tick box, `km/s` range, `±` margin, style button) | For each ticked spacecraft, and **Earth / L1** (Wind, ACE), a region of the Sun it can be magnetically connected to: its latitude, and the longitudes a Parker spiral for the speed range (default 300-800 km/s, line with the speeds at its ends) reaches at 2.5 Rsun, **west** of the sub-spacecraft point, widened by the margin (default ±10°); drawn radially on the disk, far-side parts at the limb. The panel lists the longitudes (Stonyhurst and Carrington). No PFSS: in a coronal hole the photospheric footpoint can lie elsewhere. Typical shifts: PSP near perihelion 1-3°, Solar Orbiter at 0.3 au 9-24°, L1 30-80°. The **style** button cycles: *outline* (default: coloured border with the spacecraft's name), *shaded* (translucent fill), *dim outside* (outlines, and the image outside all regions darkened). |
| save png / save pdf | Save the image as you see it (zoom, contrast, stretch, marked points), without the widgets, next to the frame's .fits file with the same name: `plume_002/171/plume_002_171_20190502_180009.png`, ... (overwrites). |
| toolbar | Zoom / pan as usual; clicks are ignored while zoom or pan is on. |
| `< prev` / `next >` (right), keys `b` / `n` | Only with several frames in the folder: step through them in time order. Image, grid, points (kept in arcsec), inset (same place) and spacecraft views follow; contrast, stretch and zoom are kept. The next frame (in the stepping direction) is read in the background; a step takes ~2 s. Clicks or keys pressed while a step is loading are ignored (the window title shows "loading frame ..."), so they do not pile up. |
| close the window | The **last** point is saved to `ref_coords.csv` (replacing that plume's earlier row, if any), with the DATE-OBS of the frame it was clicked on. Nothing is saved if no point was clicked. |

The image is the full disk, exposure-normalized (DN/s), with the limb and a
Stonyhurst grid every 15°: latitudes are labelled just outside both limbs and
longitudes along the equator; when zoomed in, the top / right axes also show
Stonyhurst longitude / latitude ticks where grid lines reach the frame. Zoom in on the plume before clicking for an accurate position.

## Planet and spacecraft positions: `plot_positions.py`

```bash
./plot_positions.py                                  # now
./plot_positions.py --date "2022-12-11 13:16" --step 1d
```

| Control | What it does |
|---|---|
| date box | Type a date (UTC) and press Enter. |
| `<<` / `>>` | Step back / forward by the selected step. |
| step | 1h, 4h, 12h, 1d, 3d, 1wk, 1mo (= 30 days). |
| now | Go to the current time. |
| Parker spirals + wind speed (km/s) | Field lines from each shown spacecraft (and Earth) down to 2.5 Rsun for the slowest (dashed) and fastest (dotted) speed of the range (default 300-800), the band between them shaded; also in the saved figures. |
| table tick boxes | Show / hide an object on the plot. At the start: Sun, Earth, Solar Orbiter, PSP. |
| save png / save pdf | Save a figure of the ticked objects only: the plot scaled to the farthest of them (+10%), their table rows and a legend, as `positions_YYYYMMDD_HHMM.png` / `.pdf` in this folder. |

The plot is in Stonyhurst heliographic coordinates: Sun at the centre, Earth
always at longitude 0. The table gives each body's distance from the Sun (AU
and solar radii) and Stonyhurst longitude / latitude (longitude = angle from Earth, + west).

Planets come from the astropy ephemeris. Spacecraft come from JPL Horizons, in
120-day blocks downloaded one spacecraft at a time (Horizons refuses parallel
requests), and only for the ticked spacecraft: the first view of a period takes
about 4 s per spacecraft, then stepping within it is instant. Every block downloaded is added to `cache/horizons_positions.json` and read from there next time (blocks that were predictions are downloaded again after 30 days). A spacecraft shows
"not available" before its launch or where Horizons has no data.

## Using the result

`ref_coords.csv` gives what the light-curve pipeline needs for a new plume in
`~/projects/project3/light_pipeline/plumes.csv`:

| `ref_coords.csv` | `plumes.csv` |
|---|---|
| `xcen`, `ycen` | `xcen`, `ycen` |
| `date_obs` (or the `--ref-date` you used) | `t_ref` (format `YYYY-MM-DD HH:MM:SS`) |

`start_date`, `end_date`, `fovx`, `fovy` and `cadence_min` are still filled in
by hand (from the movie). Note the plume names: the movies use three digits
(`plume_002`) and the pipeline uses four (`plume_0002`).

## Environment

Both scripts run with the light_pipeline Python environment (their first line
is `#!/home/biswanath/projects/project3/light_pipeline/.venv/bin/python`) and
take the email list from `~/projects/project3/light_pipeline/mails.py`. The
system `python3` can't be used: its astropy is built for numpy 1 but numpy 2
is installed. If `light_pipeline` is moved or renamed, update the first line
of both scripts and `LIGHT_PIPELINE` in `download_ref_file.py`.

To run them from any folder as `download_ref_file` / `plot_ref_file` (output
still goes to this folder):

```bash
ln -s ~/filtered_movies/download_ref_file.py ~/.local/bin/download_ref_file
ln -s ~/filtered_movies/plot_ref_file.py ~/.local/bin/plot_ref_file
```

## Solar Orbiter: search and fields of view

No FITS files are needed to find out what Solar Orbiter observed: everything comes from
the SOAR TAP service (https://soar.esac.esa.int/soar-sl-tap/tap), queried by month and
kept in `cache/soar/`.

Products (`--products`): `hri174` (HRIEUV 174 Å), `hriopn` (HRIEUV with the filter wheel open:
same 174 Å passband, more light; e.g. the daily 2-h, 10-s cadence polar coronal hole runs of
October 2025, mostly without a SOOP name; pointing from the FITS header table, as it is not
in SOAR's FOV table), `hrilya`, `phi-hrt`, `spice` (pointed: FOVs), `fsi174`, `fsi304`,
`phi-fdt` (full disk).

### The cache

Everything queried or downloaded goes into `cache/` and is read from there next time:

| Folder / file | What | Queried again |
|---|---|---|
| `cache/soar/` | SOAR metadata, one `.csv.gz` per kind and month (`items_*`, `fov_*`, `hdr_*`, `spicemode_*` = SPICE study of each file); all months since 2020-02 with `./soar_meta.py --all` (also `dkist_meta.py --all`, `iris_meta.py --all`) | only with `--refresh` (`find_obs.py`), e.g. for data released since (the file date says when it was queried); incomplete months are not kept |
| `cache/dkist/` | DKIST Data Center datasets, one file per month | only with `--refresh` (`dkist_meta.py`) |
| `cache/iris/` | IRIS observations from LMSAL HCR, the full answer per month (`hcr_YYYY-MM.json.gz`; `probe/` = first test queries) | only with `--refresh` (`iris_meta.py`), e.g. for recent days whose level-2 entries were not ready |
| `cache/soops/soops_summary.csv` | the SOOP summary table | when `soops-summary.html` is replaced by a newer copy |
| `cache/soops/pages/` | SOOP description pages (aim, instrument modes), one JSON per SOOP; all 27 that ESA has were fetched on 2026-11-04 (`./soops.py --pages`; the calibration and synoptic SOOPs have no page) | never (`./soops.py --pages --refresh`) |
| `cache/spice/studies.json` | SPICE studies: description, raster, spectral windows | never (`./spice_meta.py ... --refresh`) |
| `cache/horizons_positions.json` | spacecraft trajectories (JPL Horizons), 120-day blocks | never, except blocks that were still predictions when downloaded (after 30 days) |
| `cache/aia_synoptic/` | AIA synoptic images (1024², ~1 MB) used by `plot_fov.py` | never |
| `cache/fsi/` | EUI/FSI L2 images (~4 MB) used by `plot_fov.py` and `browse.py` | never |
| `cache/web/` | AIA / FSI images reduced for `browse_web.py`: `<file>_1024.png` (north up, default display) + `.json` (extent in arcsec, observer, colour map, unit) + `.npy` (the reduced values, float16, for other displays) | never (delete to remake) |
| `cache/search/v1/` | Keyword-search index (`search_meta.py`), one file per source and month | rebuilt for a month when its metadata in the cache is newer |
| `cache/fov_search/` | Common-FOV index (`fov_search.py`): time, observer, FOV corners per observation, one file per source and month | rebuilt for a month when its metadata in the cache is newer |
| `cache/aia_full/<wavelength>/` | full-resolution AIA frames ("original data" in `browse.py`) | never |
| `cache/solo/<product>/` | HRIEUV / HRI-Lya / PHI-HRT files ("original data" in `browse.py`) | never |

EUI files are downloaded from the EUI team's server at ROB (https://data.observatory.be/science/solo/eui/,
the same files, about 6x faster than SOAR), other files from SOAR.

Deleting a file or folder there only means it is fetched again when needed.

**All metadata is cached** (bulk download 2026-11-04): SOAR every month since 2020-02 (119 MB),
DKIST since 2020-01 (4071 datasets; refreshed 2026-11-04 to add the experiment descriptions and target types), IRIS since 2013-07 (87 MB, ~55,700 observations).
So any date opens without a metadata query. Known gaps: IRIS has no observations 2018-10-27 to
mid-December 2018 (checked in HCR, not a download failure); Solar Orbiter has no remote-sensing
files in parts of the 2020-21 cruise phase.

**Updating** (the cache is never refreshed by itself; recent months were saved as they were:
on 2026-11-04 SOAR 2026-05 to 2026-11 were not yet released, IRIS and DKIST 2026-11 partial):

```bash
./soar_meta.py --from 2026-05 --to 2026-11 --refresh   # Solar Orbiter: all browser products, FSI, SPICE studies
./iris_meta.py 2026-10-01 2026-11-30 --refresh         # IRIS (level-2 entries of recent days come later)
./dkist_meta.py 2026-11-01 2026-11-30 --refresh        # DKIST
./soar_meta.py --all; ./iris_meta.py --all; ./dkist_meta.py --all   # months added since (cached ones skipped)
```

Then restart `browse_web.py` (it keeps loaded months in memory). `find_obs.py --refresh`
refreshes only the products it is given (default HRIEUV and PHI-HRT), not FSI or the SPICE
studies, so use `soar_meta.py` for the browsers.

Beware: long SOAR TAP responses are often cut off without any error (on 2025-11-03 a
month of HRIEUV rows, 34k, came back with 21k-29k rows in six tries; the rows come in
no particular order, so the loss is scattered over the month, not at its end). So every
query is checked against a `COUNT(*)` with the same selection and repeated if short;
months with more than 4000 rows are queried day by day; an incomplete month is not
cached (a warning says so: run again). The first query of a busy month takes a few
minutes; after that it is read from the cache.

```bash
./find_obs.py 2025-10-01 2025-10-31                        # HRIEUV + PHI-HRT blocks
./find_obs.py 2025-01-01 2026-03-31 --products phi-hrt --soop Polar
./find_obs.py 2025-01-01 2025-12-31 --min-lat 60           # FOV reaches |lat| >= 60
./find_obs.py 2025-01-01 2025-12-31 --both hri174 phi-hrt  # HRIEUV and PHI-HRT together
./find_obs.py 2025-03-01 2025-04-30 --products hri174 phi-hrt fsi174 spice --csv blocks.csv

./plot_fov.py 2025-10-20T12:30                 # FOVs on AIA 193 (Earth view) + FSI 174 (SolO view)
./plot_fov.py 2025-03-16T12:00 --aia 171 --fsi 304 --zoom
./plot_fov.py 2025-10-20T12:30 --positions      # + panel with the spacecraft positions and Parker spirals
./plot_fov.py 2025-10-30T19:30 --window 24h --connect PSP "Solar Orbiter" "Earth / L1" --speeds 350-750
./dkist_meta.py 2025-10-01 2025-10-31          # what DKIST observed (bounding boxes, lines)
./iris_meta.py 2025-10-10 2025-10-17 --pole    # what IRIS observed at the poles
```

Products (`--products`): `hri174` (HRIEUV, normal and short exposures), `hrilya`,
`phi-hrt` (one row per PHI-HRT inversion, the `blos` file; `bazi`, `binc`, `icnt`,
`vlos`, ... come with it), `fsi174`, `fsi304`, `phi-fdt`, `spice` (times only).

A **block** is a run of files of one product with the same SOOP, no gap longer than 3 h
and no re-pointing by more than half the FOV. The `find_obs.py` table gives, per block:
start, duration, number of files, median cadence, SOOP; Solar Orbiter's Stonyhurst lat /
lon (lon = angle from Earth) and distance; the Stonyhurst lat / lon of the FOV centre,
the fraction of the FOV on the disk and the latitudes it covers; the fraction of that
on-disk part Earth also sees and mu seen from Earth at the centre (how close to Earth's
limb it is); and where the pointing came from.

Pointing (in `soar_meta.py`): SOAR has two FOV tables, `v_eui_hri_fov` and
`v_phi_hrt_fov`, giving two opposite corners of each image (as seen from Solar Orbiter)
and its roll angle, from which the four corners follow; and the FITS-header tables
`v_eui_sc_fits` / `v_phi_sc_fits` (CRVAL, CRPIX, CDELT, PC, NAXIS). The two agree to
better than 0.1" (checked on 14 386 HRIEUV and 469 PHI-HRT files of 2025-03). Neither is
complete: the FOV tables miss some files (56% of the HRIEUV files of 2025-10), the PHI
header table misses whole months (2024-10, 2025-09, 2025-10). So: FOV table first, then
the headers (slow, ~1 min per month for EUI, only for months with files still missing),
then the nearest file of the same block (`inferred` in the `src` column).

`plot_fov.py`: AIA background from the JSOC synoptic archive (1024², level 1.5, 2-min
cadence, ~1 MB, plain HTTP: no export to wait for); FSI from SOAR (~4 MB; the nearest
file within 1 h, else 6 h, else 1 day: the time used is in the panel title).
- Solar Orbiter FOVs: on FSI as observed (dashed off the limb). On AIA, the part on
  Solar Orbiter's disk as its **footprint on the solar surface** (solid where AIA sees
  it, dotted behind the disk). The part beyond Solar Orbiter's limb sees the corona: for
  HRIEUV, HRI-Lya and SPICE it is drawn **dashed, going up off the limb**, on a sphere
  through the Sun's centre seen from Solar Orbiter (`plot_fov.project`; a guide, since the
  depth of the emission along the line of sight is unknown). PHI-HRT sees no corona, so
  its outline follows Solar Orbiter's limb, which near the poles runs behind Earth's limb.
  The same in `browse.py` and `browse_web.py`.
- DKIST FOVs (`dkist_meta.py`, datasets within `--dkist-window`, default 3 h): as
  observed on AIA; on FSI as the coronal Solar Orbiter FOVs on AIA (footprint, off-limb
  part dashed going up).
- IRIS FOVs (`browse_web.py` only; `iris_meta.corners`): the slit-jaw area (pale; the
  slit-jaw field moves with the slit, so it covers about raster width + one image) and the
  raster area (coral), centred on the pointing at the start, rotated by the roll; drawn
  like DKIST (as observed on AIA, off-limb part going up on FSI).
- Spacecraft (`--spacecraft`, default PSP; JPL Horizons positions shared with
  `plot_positions.py`): sub-spacecraft point and limb on both images. Solar Orbiter's
  limb is always on the AIA image, Earth's on the FSI image.
- Connectivity / source regions (`--connect`, default PSP; also `"Solar Orbiter"`,
  `"Earth / L1"`, ...): as in `plot_ref_file.py`, the longitudes a Parker spiral for
  `--speeds` (default 300-800 km/s) reaches at 2.5 Rsun, ±`--margin` (10°), at the
  spacecraft's latitude; ballistic, no PFSS. The footer lists them in Stonyhurst and
  Carrington longitude.

SOOP descriptions: https://www.cosmos.esa.int/web/solar-orbiter/soops-summary
(saved in `soops-summary.html`, read by `soops.py`; `find_obs.py` prints the runs behind its
blocks: target, coordinated observations, notes).
