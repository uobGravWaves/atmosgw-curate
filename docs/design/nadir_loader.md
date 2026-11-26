# Design: `load_nadir_sounders` (port of [`prep_airs_3d`](https://github.com/corwin365/MatlabFunctions/blob/master/DatasetSpecific/AIRS/prep_airs_3d.m))

**Status:** design of 2026-11-25, implemented in 0.4-0.5 (format version 1). Parts
still planned, not implemented: the `analysis` steps (detrend, smooth, fill)
and `data.to_netcdf` with packing.
Decisions marked *agreed* were agreed with C. Wright on 2026-11-25; anything
else is a proposal.

Companion document: [`data_formats.md`](data_formats.md), which defines the
output format this loader produces.

## 1. Purpose and scope

Port the MATLAB `prep_airs_3d` (corwin365/MatlabFunctions,
`DatasetSpecific/AIRS/prep_airs_3d.m`, as fixed in November 2026) to Python,
as `atmosgw.curate.obs.load_nadir_sounders`: load nadir-sounder swaths and
return them in a common format, on a regular grid by default.

- **Instruments now:** AIRS and IASI.
- **Instruments later**, without changing the core: CrIS, IASI-NG and similar
  cross-track scanners.
- **Products:** 3-D (retrieved temperature on levels) and 2-D (a single layer,
  e.g. brightness temperature in one or more channels). The core handles both
  now; 2-D readers come later.
- **Not in the loader:** detrending and smoothing, which become general
  analysis functions (§7), as the Curate pipeline has a separate Detrender.

## 2. Place in the Curate pipeline

From the Curate plan (*Seeing the Unseen* / GW-LIFT slides, "GW property
retrieval plan"):

```
instrument loader(s) -> Detrender -> Stockwell transform -> Parameter Engine -> Fuse
        ^                                                     ^        ^
        | characterisation package (CP)                       | CP     | geophysical context (GC)
```

What this asks of a loader, and how this design answers it:

| Plan requirement | In this design |
|---|---|
| Common interfaces and data structures ("very high" priority) | Formats specified and checked: `data_formats.md`, `data.validate` |
| Consistent treatment across instruments; extensible to new ones | One core shared with the limb loader; per-instrument readers behind a fixed contract (§5) |
| Fully incorporates instrument QC | Extensible per-variable state flags (`<variable>_flag`), shared registry of states; instruments' own flags passed through |
| Uncertainties: interfaces now, implementation later | `<variable>_uncertainty` naming; static uncertainty in the CP, linked by attribute |
| Provenance / independence for Fuse (e.g. 2-D and 3-D AIRS share radiances) | In the CP, linked by attribute |
| Sampler and CPs need sampling and viewing geometry | Native geolocation kept on request; slots for satellite geometry (§8) |
| Detrender consistent across 1-D/2-D/3-D | Detrending outside the loader, general over dimensions (§7) |

## 3. Decisions

| # | Question | Decision |
|---|---|---|
| 1 | Name | `load_nadir_sounders` (*agreed*) |
| 2 | Granules | Primary product is one continuous swath over the time range; granules are accounting. Data can be requested by granule, and every row records its granule (*agreed*) |
| 3 | Detrending, smoothing | Outside the loader, in `analysis/` (*agreed*) |
| 4 | 2-D products | No `level` dimension; a `channel` dimension where there are several (*agreed*) |
| 5 | Quality control | Extensible flags (*agreed*): one CF bitmask with a shared registry of bits (`data_formats.md` §6) |
| 6 | Uncertainty and provenance | In the characterisation package (*agreed*); data carry only per-point values supplied by the files, and a link to the CP |
| 7 | Shared helpers | Moved out of `obs/_limb_tools.py` into `obs/_common.py`, limb outputs bit-identical (*agreed*) |
| 8 | Selecting data | By `time_range` and `time_handling`, as limb; `granules=` as an option (*agreed*) |
| 9 | Formats | Written down and checked (`data_formats.md`) (*agreed*) |
| 10 | Bad points | NaN in every reader, never refilled by the loader; gaps are filled, if wanted, by a separate function before spectral analysis (*agreed*) |
| 11 | Absent quantities | "Never available" (e.g. geometry not in the files) is recorded separately from "missing here" and "bad" (*agreed*; `data_formats.md` §6) |
| 12 | Smoothing | Our own boxcar implementation, not a port of `smoothn` (*agreed*) |
| 13 | Leap seconds | Ignored: well below wave timescales (*agreed*) |

## 4. Architecture

| Layer | Limb (exists) | Nadir (new) |
|---|---|---|
| Public loader | `obs.load_limb_sounders` | `obs.load_nadir_sounders` |
| Step-by-step core (MATLAB names, compared with MATLAB) | `obs/limb.py: _load_limb_arrays` | `obs/nadir.py: _load_nadir_arrays` |
| Instrument readers | `obs/_limb_<inst>.py` | `obs/_nadir_airs.py`, `obs/_nadir_iasi.py`; later `_nadir_cris.py`, `_nadir_iasing.py` |
| Load-time helpers | `obs/_limb_tools.py` | `obs/_common.py`, shared by both |
| Array maths | `profiles.py` | `profiles.py` (vertical interpolation, unchanged) + new `swath.py` |
| Formats and flags | (implicit) | `data/formats.py`, `data/flags.py`, `data.validate` |
| Analysis | — | `analysis/detrend.py`, `analysis/smooth.py` |

### Shared with the limb loader

- **Used unchanged:** `paths.dataset_dir` and `datasets.toml` (new AIRS and
  IASI entries); `readers.find_files`, `readers.read_netcdf`; `matlab_time`;
  `geodesy` (`haversine`, `azimuth`, `reckon`, `km2deg`, which replace the
  Mapping Toolbox calls `prep_airs_3d` makes); `matlab_compat.interp1`;
  `profiles.interpolate_profiles` for the vertical grid (each swath column is
  a profile: reshape, interpolate, reshape back).
- **Moved to `obs/_common.py`:** `time_handling` rules; the load record
  (`files_read`, `files_not_found`, `multiple_matches`, variable sources and
  notes); `choose_file`; checks of shared options; the outlier rules; building
  the output Dataset and its attributes.
- **New, general:** `swath.py` (along/cross-track distance coordinates,
  regridding to a regular distance grid with the dateline unwrap, the regular
  vertical grid); solar zenith angle (port of `which_airs_retrieval`, in
  `geodesy.py` or a new `solar.py`).

## 5. Reader contract

Every nadir reader is a module with one function:

```python
load(settings, inst_info, var_names) -> (data, file_list)
```

- `settings`: the checked options, including the expanded time range and any
  requested granules.
- `inst_info`: where the files are, and the load record (as limb).
- `var_names`: the variables wanted, with MATLAB-style internal names
  (`Lat`, `Lon`, `Time`, `Temp`, `Alt`, ...) plus any `additional_vars`.
- `data`: a dict of NumPy arrays for all requested granules, in time order,
  already joined along track:
  - per point: `[cross, along]` (`Lat`, `Lon`, `Time` as datenums);
  - per point and level: `[cross, along, level]` (`Temp`, ...); 2-D products
    give `level` length 1 inside the core;
  - per level: `[level]` (`Alt`, km);
  - per row: `[along]` (`Granule`: granule number; `Day`: its day;
    `SourceFile`: 1-based index into `file_list`; `SourceRow`: 1-based row in
    that file);
  - optional geometry (§8) and instrument QC variables, same shapes.
- `file_list`: file names read, in order.

Gaps (missing granules or scans) are not closed up: the reader returns rows
in time order, and the core splits the swath into `segment`s at gaps (§6).

### Instrument registry

As `INSTRUMENTS` in `limb.py`, one entry per instrument:

| Field | AIRS | IASI |
|---|---|---|
| `reader` | `_nadir_airs` | `_nadir_iasi` |
| `n_cross` | 90 | 60 |
| `granule_seconds` | 360 | 360 |
| `granules_per_day` | 240 | 240 |
| granule origin | from the files (granule numbers) | pseudo-granules: [6(k-1), 6k) min after 00:00 UTC |
| default `n_along` per granule | 135 | 90 (rows in a pseudo-granule) |
| `time_range` | 2002-08-30 onwards (UNVERIFIED: check first file) | 2007-05-29 onwards on eepc-0184 |
| `products` | `"3d"` (Hoffmann and Alexander retrieval); `"2d"` later | `"3d"`; `"2d"` later |
| CP id (placeholder, *agreed*) | `AIRS-3D-placeholder/0` | `IASI-3D-placeholder/0` |

Adding an instrument means one registry entry, one reader and one
`datasets.toml` entry. The only constraint on a reader is that it delivers a
rectangular cross x along grid (CrIS's 3 x 3 fields of view per field of
regard must be laid out as rows and columns by its reader).

## 6. Public interface

```python
load_nadir_sounders(time_range, instrument, *,
    product="3d",
    time_handling="fuzzy",
    granules=None,
    lat_range=(-90, 90), lon_range=(-180, 180),
    original_grid=False,
    n_cross=None, n_along=None, dx_cross=None, dx_along=None,
    original_z=False, height_scale=None, height_range=(0, np.inf),
    interp_type="linear", extrapolate=False,
    keep_outliers=False,
    day_night=False,
    keep_native_geolocation=False,
    additional_vars=(), file_source=False, file_path=None,
    date_warning=True, verbose=False,
) -> xarray.Dataset
```

Options shared with `load_limb_sounders` have the same names, meanings and
defaults (`time_handling`, `lat_range`, `lon_range`, `original_z`,
`height_scale`, `height_range`, `interp_type`, `extrapolate`,
`keep_outliers`, `additional_vars`, `file_source`, `file_path`,
`date_warning`, `verbose`).

| Option | Meaning | MATLAB |
|---|---|---|
| `time_range`, `time_handling` | Which data: whole granules whose mean time is in range (as limb selects whole profiles) | `DateNum`, `GranuleId` |
| `granules` | Explicit granules: a list of `(day, number)`; overrides `time_range` selection | `DateNum`, `GranuleId` |
| `product` | `"3d"` or `"2d"` | — (3-D only) |
| `lat_range`, `lon_range` | Keep rows of the swath that enter the box (whole rows, so the grid stays rectangular) | — |
| `original_grid` | Return the native geometry, no horizontal regridding | `LoadOnly` (also skips vertical) |
| `n_cross`, `n_along` / `dx_cross`, `dx_along` | Regular grid: number of points, or spacing in km (spacing overrides number) | `NXT`, `NAT`, `dXT`, `dAT` |
| `original_z`, `height_scale`, `height_range` | Vertical grid; default: regular 3 km grid from the lowest native level, as MATLAB (its comment says "from the middle outwards", but the code starts at the bottom) | `VerticallyInterpolate`, `VerticalSpacing` |
| `interp_type`, `extrapolate` | Horizontal and vertical interpolation | `Interpolant`, `Extrapolant` |
| `keep_outliers` | Do not set flagged points to NaN (flags are always returned) | — |
| `day_night` | Add `day_retrieval` (solar zenith angle < 96 deg, the angle at which L. Hoffmann's retrieval switches between its day and night set-ups, per C. Wright) | `DayNightFlag` |
| `keep_native_geolocation` | Also return native latitude, longitude and time before regridding (for the sampler and CPs) | `KeepOldTime` (time only) |
| `additional_vars` | Extra file variables | — |
| `file_source` | Also return `source_file`, `source_row` per row and `original_files` | `Source` |

Not ported: `InputStruct` (use `swath.regrid` on your own arrays instead),
`NoISCheck`, `Python`, `PreSmooth` and `NoDetrend`/`DetrendMethod`
(analysis functions, §7), `IasiFracBad`/`IasiRegionFracBad` (replaced by the
flags and their counts in the attributes).

### Output

A Dataset in the **swath** format of `data_formats.md`, version 1:

- dimensions `(along_track, cross_track, level)`; 2-D products
  `(along_track, cross_track)` or `(along_track, cross_track, channel)`;
- `granule`, `granule_day` and `segment` as coordinates along `along_track`,
  so any row can be traced to its granule, and data can be split by granule
  for collaborators;
- data in float32 (a whole day of AIRS on the default grid is about
  90 x 32 400 x 31 values: 0.36 GB per 3-D variable in float32). Long ranges
  are best loaded a day at a time; lazy loading (dask) can come later if
  needed.

**Segments:** the swath is regridded segment by segment. A gap of more than
1.5 nominal scan intervals ends a segment; rows between segments are not
interpolated across, and the output has no rows for the gap (each segment is
regridded on its own distance axis, and segments are joined along track).
This reproduces MATLAB's joined multi-granule swath where the granules are
contiguous.

## 7. Analysis functions (outside the loader)

- `analysis.detrend.polynomial(ds, var="temperature", dim="cross_track",
  order=4)`: fits and removes a polynomial along one dimension for each line
  of the others; adds `temperature_perturbation` and
  `temperature_background`. Port of `airs_4dp_detrend` (MATLAB `Tp`, `BG`,
  Method 2; Method 1 gives the same result). Works for any dimension, so the
  same function serves 1-D, 2-D and 3-D data.
- `analysis.smooth.boxcar(ds, var, size={...})`: N-D boxcar with edge
  replication, our own implementation (the MATLAB uses `smoothn`, File
  Exchange 725, which is not ported). Its output is compared with MATLAB's
  `PreSmooth` results as a check.
- `analysis.fill.fill_missing(ds, var, dims=..., method="linear")`: fills NaN
  points from their neighbours, for spectral methods that cannot take gaps,
  and sets the `filled` state on every point it fills (`data_formats.md` §6).
  The loader never calls it: filling is a visible, separate step.

## 8. Geometry hooks

The swath format has optional variables for viewing geometry
(`satellite_latitude`, `satellite_longitude`, `satellite_altitude`,
`view_zenith_angle`, `view_azimuth_angle`, `scan_position`), each
`[cross, along]`. Readers fill those the files hold; the regridding
interpolates them like latitude and longitude. A variable the files do not
hold is not in the output, and the dataset's `availability` attribute says
why (`not_in_files` here), so it can never be mistaken for a bad or missing
value (`data_formats.md` §6). The AIRS files on eepc-0184 had `l1_sat_*`
removed to save space; re-fetching the raw files will let the AIRS reader
fill them.
`scan_position` (the native cross-track index) is always available and is
returned with `keep_native_geolocation`.

## 9. Processing order

As MATLAB, so the step-by-step core can be compared with it:

1. check options; expand the time range; find granules;
2. read (reader): AIRS daily or granule files; IASI orbit files overlapping
   each pseudo-granule; time from seconds since 2000-01-01 (leap seconds
   ignored, as MATLAB: agreed, well below wave timescales); file fill values
   become NaN with the `missing` state; so do the undeclared fill values and
   bad times of the AIRS and IASI files (`obs._common.nadir_fill_to_nan`:
   -999, -9999, 0/0 positions, zero, stray and misdated times), whose
   positions and times are then `bad_geolocation`; scans with any such
   position or time are their own segments, not regridded;
3. IASI screen (`obs._screen`; default rolling 5 x IQR) and the 100-400 K range check, setting
   their states; flagged points become NaN. (MATLAB refills IASI points from
   the nearest along-track value; Python does not, by design. The MATLAB
   comparison is therefore made on the screened values before MATLAB's
   refill, and separately on `fill_missing` output with
   `method="nearest"` along track.)
4. missing heights: the AIRS fallback height grid, as MATLAB (flagged in a
   variable note);
5. split into segments; regrid each horizontally (distance space, dateline
   unwrap) unless `original_grid`;
6. vertical grid unless `original_z` (5 km native levels above 60 km are
   interpolated to 3 km, as MATLAB);
7. outliers to NaN unless `keep_outliers`; geographic selection;
8. day/night, file source, attributes; validate against the swath format.

## 10. Validation

- **MATLAB reference:** `test_prep_airs_3d.m` (corwin365/MatlabFunctions-Tests)
  and its reference outputs (`prep_airs_3d_reference.mat`, made from
  MatlabFunctions commit 059037e; to be remade from the commit with the
  November 2026 fixes once merged). The private layer `_load_nadir_arrays`
  is compared case by case, as the limb port. Where Python differs by design
  (continuous swath across gaps, float32 output), the comparison is on the
  matching cases (contiguous granules, float64 inside the core).
- **Known answers:** the MATLAB test's checks, ported (midnight crossing, day
  shift, dateline, longitude rotation, IASI bins and file starts, outlier
  screen on a synthetic file).
- **Format:** every output passes `data.validate(ds, "swath")`.
- Tolerances stated in each test, as `CONTRIBUTING.md` requires.

## 11. Checkpoints

Work happens in a separate worktree (`~/Code/atmosgw-curate-nadir`), on its
own branches, so it can run alongside the limb-loader work.

0. This design and `data_formats.md` (this pull request).
1. `data/` (format spec, flag registry, validator) and `swath.py`, solar
   zenith angle; known-answer tests. *New files only.*
2. AIRS 3-D reader and `load_nadir_sounders`, matching the MATLAB reference.
3. *Deferred (agreed 2026-11-25), waiting for a general Detrender:*
   `analysis.detrend`, `analysis.smooth` (MATLAB `Tp`, `BG`). Until then the
   loader returns temperatures only, with no perturbations.
4. IASI reader. *Done.* The outlier screen is swappable (`obs._screen`;
   `load_nadir_sounders(screen=...)`): IASI uses a rolling 10-min window,
   day and night separately, 5 x IQR (agreed 2026-11-25, from measurements);
   AIRS has none (its 3 x IQR flags were real variability); MATLAB's
   pooled screen is kept for comparison. Regridding defaults to 90 x 135 per
   granule, as MATLAB (native 60 x 90).
5. The 2-D path, `data.to_profiles`, docs. *Done* (release left for C. Wright).
   2-D products work in the core (tested with a stand-in reader, and on
   L. Hoffmann's AIRS 2-D perturbation files in a throwaway reader). Those
   files (`airs_YYYY_DDD_GGG.nc`, 135 x 90 per granule; e.g. 2020 days 20-60
   on the group's Linux machine) hold
   brightness temperature in channel groups `bt_4mu`, `bt_8mu`,
   `bt_15mu_low`, `bt_15mu_high`, each with a perturbation (`_pt`) and a
   variance (`_var`): an AIRS 2-D reader must decide which to load, the
   channel wavenumbers, and how to label the file's own perturbations
   (already detrended by their producer).
   *AIRS 2-D reader done* (`obs/_nadir_airs2d.py`), with these choices
   (proposals): all four channel groups by default (`channels=` to choose);
   nominal wavenumbers 10^4 / wavelength (UNVERIFIED) and a `channel_name`
   label, as both 15 micron groups share one; the files' perturbations and
   variances returned as `brightness_temperature_perturbation` and
   `brightness_temperature_variance`, `origin` "derived" (NaN for 8.1 micron,
   which has none); files in `YYYY/DDD/` or `DDD/`; data folder `AIRS` `2d`
   in `datasets.toml`.
6. Move the shared helpers into `obs/_common.py` and use them in both loaders
   (limb outputs bit-identical). *Done* (2026-11-25): `obs/_common.py` holds the
   time rules, the load record and file choice, the shared option checks,
   the temperature range check and the load attributes; limb and nadir
   outputs identical on real data (32 cases).

## 12. Open items

- Start dates of AIRS and IASI data, for the registry. UNVERIFIED.
- **AIRS daily files store time as float32**: `l1_time` (seconds since 2000)
  comes in steps of 16 s in 2007 and 64 s in 2025, so times are up to 8-32 s
  off (checked 2025-10-08, granule 93: up to 32.0 s against the granule
  file). The swath is therefore split into segments only at gaps over 180 s.
  Like the temperature rounding, rebuilding the daily files would fix it.
- **AIRS daily files are rounded to whole kelvin.** `merge_airs_into_daily.sh`
  stores `ret_temp` with `short()` and no scale factor. Checked on one
  profile (2025-10-08, granule 93): granule file 279.215, 267.145, 248.637 K;
  daily file 279, 267, 249 K (rounded to nearest, so unbiased, but with
  about 0.29 K rms quantisation noise). The original granule files keep full
  precision. *Agreed:* keep the option of replacing the daily files with
  unrounded ones (about 0.01 K, e.g. int16 with `scale_factor` 0.01, or
  float32). The reader must not assume whole kelvin: `readers.read_netcdf`
  already applies `scale_factor` and `add_offset`, so either form reads the
  same way, and nothing downstream may depend on the 1 K steps. The
  characterisation package records the precision of the files in use.
- *Agreed:* applying the format spec to `load_limb_sounders` output
  (attributes, state flags and `availability`; values unchanged) is a separate
  change, after the
  limb work in progress is merged. *Done* (2026-11-25).

## 13. MATLAB -> Python names

`l1_lat` latitude, `l1_lon` longitude, `l1_time` time, `ret_temp`
temperature, `ret_z` altitude, `ret_press` (IASI) via `additional_vars`,
`Tp` temperature_perturbation, `BG` temperature_background, `DayNightFlag`
day_retrieval, `Spacing` the `grid_spacing_km` attribute, `Source` /
`FilePath` original_files.
