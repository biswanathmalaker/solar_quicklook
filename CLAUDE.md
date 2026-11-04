# solar_quicklook

Quick-look and search tools for coronal and Solar Orbiter data: SDO/AIA, Solar Orbiter (EUI, PHI, SPICE) from SOAR, DKIST, SOOPs and spacecraft positions. They serve the polar-plume paper in `~/projects/project5` (see its CLAUDE.md and PLAN.md). The folder was copied from `~/quick_aia` on 2026-11-03. Leave `~/quick_aia` and `~/solar_orbiter` untouched.

## Rules from the user
- **Python:** `~/projects/project3/light_pipeline/.venv/bin/python` (sunpy 6, numpy 2). Scripts carry it as their shebang. Don't create another environment; the system python3 doesn't work.
- **Writing and running code:** allowed freely here, as long as it is research-relevant. Dig for details the user might miss.
- **Portable:** every path is relative to this folder, and nothing is stored outside it. Copy the folder (code + `cache/`) and it works offline for everything already fetched (tested 2026-11-04 with the network blocked). The only absolute path is the venv in the shebang lines.
- **Cache:** everything queried or downloaded goes in `cache/` and is never fetched again automatically. Re-querying happens only with an explicit `--refresh`. Incomplete results (e.g. a truncated TAP month) are not cached.
- **Full metadata in the cache (user's choice; done 2026-11-04, no failures):** `./soar_meta.py --all` (since 2020-02, 119 MB), `./dkist_meta.py --all` (since 2020-01, 3698 datasets), `./iris_meta.py --all` (since 2013-07, ~55,700 observations, 87 MB) cache every month not yet cached, so browsing any date needs no metadata query.
  - Recent months were cached as they were: SOAR 2026-05..11 not yet released (Sep-Nov empty), IRIS / DKIST 2026-11 partial. Update with `./soar_meta.py --from 2026-05 --to 2026-11 --refresh`, `./iris_meta.py START END --refresh`, `./dkist_meta.py START END --refresh`, then restart `browse_web.py` (months are kept in memory). `find_obs.py --refresh` only covers the products given to it.
  - Real gaps: IRIS none from 2018-10-27 to mid-Dec 2018 (checked with smaller HCR queries); Solar Orbiter none in parts of the 2020-21 cruise phase.
- **commands.md:** one line per command or option, nothing extra. Keep it, README.md and `help.py` in step with the code.
- **Colours:** every instrument's field of view has its own colour (asserted in `plot_fov.py`). The two HRIEUV products share lime.
- **Summaries:** say in words what each instrument measures, with time range and cadence. Bare numbers are not enough (user feedback).

## Files
| File | Role |
|---|---|
| `soar_meta.py` | SOAR TAP metadata, cached by month in `cache/soar/`. Main functions: `files()`, `blocks()` (observing blocks), pointing, `spice_modes()`, `download()` (EUI from the ROB mirror first). |
| `find_obs.py` | Command-line table of observing blocks, with filters for SOOP, latitude and instrument overlap. |
| `plot_fov.py` | Static figure: fields of view on AIA (surface footprints, dotted behind the limb) and FSI; PSP and Parker-spiral source regions; DKIST. |
| `browse_web.py` | **The web browser** (Panel 1.9 + Bokeh 3.9 on localhost:5006/browse, built 2026-11-03).<br>Timeline of blocks, AIA and FSI with zoom and hover, tabs (summary, files, DKIST).<br>Reuses `browse.MetaStore`, `browse.summary_lines` and the other helpers, so don't delete `browse.py` yet.<br>Images are reduced once into `cache/web/`; FOVs use closed-form `hgs_to_hpc` (checked against sunpy to 1e-10″).<br>A step takes ~0.1 s once cached, 5–7 s cold. |
| `browse.py` | Old matplotlib browser (12–20 s per step). Still holds the shared step and summary code. |
| `iris_meta.py` | IRIS observations from LMSAL HCR, the raw answer cached per month in `cache/iris/hcr_YYYY-MM.json.gz` (~40 s and ~0.5 MB a month, once). `observations(t0, t1)` gives one row per observation; `corners(r, "raster" / "sji")`; `raster_text`, `sji_text`. |
| `dkist_meta.py` | DKIST Data Center API (no login); cached in `cache/dkist/`. |
| `soops.py` | ESA SOOP summary (`soops-summary.html`, saved by hand) parsed into `cache/soops/`. `page(soop)` gets a SOOP's aim and instrument modes; cached in `cache/soops/pages/`. |
| `spice_meta.py` | SPICE study description and spectral windows, read from FITS headers by HTTP byte range on the IAS release server; cached in `cache/spice/studies.json`. |
| `plot_positions.py` | Spacecraft positions (JPL Horizons, cached in `cache/horizons_positions.json`). |
| `download_ref_file.py` / `plot_ref_file.py` | Download full-resolution AIA frames (JSOC; `folder=` argument) and click plume coordinates. |

## Non-obvious facts about the data
- **SOAR TAP:** https://soar.esac.esa.int/soar-sl-tap/tap/sync, ADQL, CSV.
  - Long responses are silently truncated. `fetch_verified` checks every query against COUNT(*); months with more than 4000 rows are queried day by day.
- **Pointing tables are incomplete:**
  - The FOV tables (`v_eui_hri_fov`, `v_phi_hrt_fov`) give the corners of the first and last pixel plus crota.
  - The header tables are `v_eui_sc_fits`, `v_phi_sc_fits` and `v_spi_sc_fits` (SPICE has no FOV table).
  - Use the FOV table first, then headers. Files with no pointing take it from the nearest file in their block ("inferred").
- **HRIEUV open filter** (`eui-hrieuvopn-image`, product `hriopn`): still 174 Å. This is where the high-cadence runs are, e.g. the daily 10-s polar coronal hole runs of Oct 2025. It isn't in the FOV table.
- **SOAR SPICE headers** give the study (`observation_mode`) but only the first spectral window. Full windows: `spice_meta.py`.
- **EUI files:** the ROB mirror `data.observatory.be/science/solo/eui/L2/YYYY/MM/DD/` is byte-identical to SOAR and about 6× faster.
- **AIA synoptic images** (1024², 2-min cadence, level 1.5) come from `jsoc1.stanford.edu/data/aia/synoptic/`.
- **IRIS:** the iris.lmsal.com/search page is a front end to the LMSAL HCR API. No login is needed.
  - `https://www.lmsal.com/hek/hcr?cmd=search-events3&outputformat=json&instrument=IRIS&startTime=…&stopTime=…` returns one record per observation. This includes the aim (`sciObjectives`), `planners`, `descriptions`, pointing and FOV, `spectral_window_names`, data URLs, and `groups` (one per raster or SJI channel: cadence, steps, exposure, image count).
  - Per observation that is ~16 kB. `outputformat=jsonMedium` is ~0.4 kB but lacks the aim and planner.
  - A month (~380 observations) is not truncated (`iris_meta` asks `limit=5000` and refuses to cache if reached).
  - FOVs: `hpc_geom_inner` = the raster box (`raster_fovx × raster_fovy` on xCen, yCen). `xFov × yFov` = the area the slit-jaw images cover over the raster (≈ raster width + one SJI width), which equals the SJI groups' fovx/fovy. All 2025-10 rolls are 0; the roll sign in `iris_meta.corners` (counter-clockwise) is unchecked.
  - In the Raster group, HCR gives sit-and-stare runs as N "steps" of 0″ (N = exposures) and signed step sizes; `basic_cadence_avg_asrun` = per step, `cadence_avg_asrun` = per raster.
- **FOVs seen from the other vantage point** (`plot_fov.project`, decided with the user 2026-11-03):
  - The on-disk part is drawn as its surface footprint: solid where seen, dotted behind the disk.
  - The off-limb part of coronal instruments (HRIEUV, HRI-Lya, SPICE; DKIST on FSI) is drawn **dashed, going up off the limb**, on a sphere through the Sun's centre, closed form (= sunpy `SphericalScreen`, checked to 1e-8″).
  - PHI-HRT (photosphere, no off-limb data) instead follows Solar Orbiter's limb ("wrap"). Don't wrap the coronal ones.
- **Observing blocks** start a new block when there is a gap of more than 3 h, the SOOP changes, or the pointing moves by more than max(60″, half the field of view's longer side).

## Status and next work (updated 2026-11-04)
**`browse_web.py` is done and the user is happy with it.** The user tests it in their own browser and reports problems. **Don't drive Chrome** to check it (the user said so); test headlessly by driving `App` and `poll()` from Python.
- **Features:**
  - timeline (incl. IRIS);
  - AIA and FSI with FOVs (off-limb rule above) and hover on FOVs;
  - summary, files, DKIST, IRIS and positions tabs;
  - buttons: positions, download original data, show Solar Orbiter originals;
  - arrow keys work (the user confirmed);
  - a cursor readout under each image: helioprojective x, y and Stonyhurst lon, lat and μ, computed in JS (`READOUT_JS`).
- **μ in the readout is exact** for a finite observer distance. `soar_meta.describe_fov` / `cos_angle` treat the observer as infinitely far, so μ comes out ~0.005–0.015 high near the pole at 0.66 au. Switch them to the exact form before quoting μ in the paper.
- **Next, in order:**
  1. **IRIS: done 2026-11-04** (`iris_meta.py`; in `browse_web.py` a timeline row, slit-jaw and raster boxes on AIA and FSI, summary lines, IRIS tab with level-2 links). Not done: downloading IRIS data (level-2 tarballs are 0.1–1.5 GB; links only).
  2. **Custom image display** (after IRIS; the user will specify): choices beyond the current default for AIA, FSI and the originals (vmin/vmax, stretch, colour map). The current default is in `light_image`: asinh stretch with a 99.7% percentile interval, the map's own colour map, and PHI B_LOS grey ±100 G. The reduced PNGs in `cache/web/` bake that default in, so custom display probably needs the reduced data array, not only the RGBA image.
  3. **Local store** at 1-h cadence for fast lookup: AIA images; Solar Orbiter summaries and metadata; DKIST fields of view and metadata only (no DKIST data, too large); all IRIS records. Before bulk downloads, confirm the date range, the AIA channels and resolution, and the format.
  4. **More buttons:** the user will specify them.
