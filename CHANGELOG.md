# Changelog

All notable changes to this package. Versions follow semantic versioning;
the package stays at 0.x until the ports match the MATLAB
reference outputs.

## [Unreleased]
- IASI: temperatures of exactly 100 K, the retrieval's lower bound (the files
  never go below it, and it fills whole bad profiles), are now
  `outside_physical_range` and removed, whatever screen is used (2010-07-01:
  30,805 points in 6,467 profiles). Limb and AIRS unchanged (MLS has a few
  values at exactly 100 K, kept as MATLAB). [Claude]

## [0.4.0] - 2026-11-25
- One set of data-state rules for both loaders and `to_profiles`
  (`docs/design/data_formats.md` section 6.3), checked by the new
  `tests/test_states_conformance.py`. Limb: a NaN between a profile's first
  and last value is now `missing`, not `no_observation`. Nadir: points
  outside a segment's own coverage get `extrapolated` (with `extrapolate`).
  Values unchanged. [Claude]
- New dataset attributes in both formats: `preprocessing` (what was done to
  the files before loading; `data` is then "preprocessed": the GNSS files,
  pre-interpolated by C. Wright, and AIRS 3-D loads from the daily files, not
  the granule files) and `heights`. Written by both loaders; not required by
  `validate` in format version 1, so older files stay valid. [Claude]
- Nadir: `source_variable`/`comment` attributes from the readers' record, as
  for limb; `state_counts` covers every data variable. [Claude]
- README: real status, both loaders, and a validation-status section;
  `docs/usage.md`: state table, "same option, different behaviour". [Claude]
- `to_profiles` keeps the file provenance of a swath loaded with
  `file_source=True` (`source_file` over profiles and levels, `source_row`
  per profile, `original_files`), lists `source_profile` as `not_provided`,
  and keeps the swath's reason for other absent variables. [Claude]

## [0.3.6] - 2026-11-25
- `load_limb_sounders`: pressure estimated from height is marked `derived`
  point by point, and `origin` ("derived", "mixed", "retrieved") goes by the
  profiles in the output, not the files read; with `extrapolate=True`, points
  outside a profile's own heights have the `extrapolated` state. Values
  unchanged. [Claude]

## [0.3.5] - 2026-11-25
- New usage guide, `docs/usage.md`: loading limb profiles and nadir swaths,
  times, formats, data states, availability and saving. Loader docstrings
  brought up to date (IASI screen, segment rule). [Claude]
- `load_limb_sounders` output is now in the common "profile" format: a state
  flag (`<variable>_flag`) and an `origin` for each data variable, and the
  attributes `format`, `format_version`, `product`, `characterisation`,
  `availability` and `state_counts`. Values are unchanged (checked on real
  data for eight instruments). The profile format accepts float64 data
  variables. [Claude]
- Shared loader helpers (time rules, load record, file choice, option checks,
  range check, load attributes) moved into `obs/_common.py` and used by both
  `load_limb_sounders` and `load_nadir_sounders`; outputs unchanged. [Claude]
- Faster nadir regridding: `swath.regrid` finds cells and weights once per axis
  and applies them to every field, level and flag (in blocks), replacing
  scipy's interpolator; vertical interpolation onto shared heights is one
  array operation; data states are carried by OR-ing the flags of the points
  used (new `flags=` in `swath.regrid`, `swath.interpolate_level_flags`)
  instead of interpolating an indicator per state. A regridded day: AIRS
  27 s -> 12 s, IASI 96 s -> 23 s and 4.6 -> 2.8 GB peak. Results identical
  except times, which differ by up to 31 microseconds (datenum rounding).
  [Claude]
- `load_nadir_sounders`: AIRS 2-D (`product="2d"`): L. Hoffmann's granule files,
  brightness temperature in four channel groups (`channels=` to choose), with
  the files' perturbations and variances as derived variables. Readers and
  data folders can now differ by product. [Claude]
- `load_nadir_sounders`: 2-D products (brightness temperature, no levels,
  one or several channels) handled by the core; no 2-D reader yet. [Claude]
- New `data.to_profiles`: a 3-D swath as profiles (profile format), with
  pressure estimated from altitude and marked derived. [Claude]
- New `data` package (proposal, format version 1): the common "swath" and
  "profile" formats, a registry of data states (why a value is absent, or
  what kind of value it is) recorded per variable as CF flags, an
  `availability` record for variables that do not exist in a Dataset, and
  `data.validate`. See `docs/design/data_formats.md`. [Claude]
- New `swath` module: distances across and along a swath, regridding onto a
  regular distance grid (with the dateline unwrap) and the regular vertical
  grid; port of the regridding in the MATLAB `prep_airs_3d`, matching it on
  synthetic swaths (within 1e-9). [Claude]
- New `solar` module: solar zenith angle and the AIRS day/night retrieval
  switch (port of `which_airs_retrieval.m`). [Claude]
- New `obs.load_nadir_sounders` (preliminary; AIRS 3-D): nadir-sounder swaths
  in the swath format, port of the MATLAB `prep_airs_3d` without its
  detrending; continuous swath with granule and segment tracking, data states,
  optional day/night, native geolocation and file source. Matches the MATLAB
  reference outputs for every AIRS case (within 3e-11 degrees, 1e-12 K,
  4e-10 days; load-only bit-identical). AIRS added to `datasets.toml`.
  [Claude]
- `obs.load_nadir_sounders`: IASI 3-D (orbit files, 6-min pseudo-granules,
  every file overlapping a granule; files without temperature skipped and
  recorded). Outlier screens are swappable (`screen=`: a name, settings or a
  function; recorded in the `screen` attribute): IASI uses a rolling 10-min
  window, day and night separately, 5 x IQR, with data padded so flags do not
  depend on the request; AIRS has none; MATLAB's pooled screen is kept.
  Screened points are NaN and flagged, never refilled. Matches the MATLAB
  reference (with the pooled screen) wherever Python has a value. IASI added
  to `datasets.toml`. [Claude]
- Design notes for the nadir-sounder loader and the data formats
  (`docs/design/`). [Claude]
- `tests/matlab/sweep/`: a by-hand check of `load_limb_sounders` against the
  MATLAB `get_limbsounders` on many dates (every instrument and year with
  data; first run 2026-11-25: 1420 calls, all identical or within stated
  last-bit tolerances). Not run by `pytest` or CI. [Claude]
- `tests/test_limb_matlab.py`: four more MATLAB reference cases found by the
  sweep (MLS heights kept, ALOMAR heights to 120 km, SABER `dx`, an MLS day
  with v05 and v06 files), with their stated tolerances; new reference file
  (the 33 earlier cases bit-identical). [Claude]

## [0.3.4] - 2026-11-25
- `load_limb_sounders` needs about a third less memory at its peak (e.g.
  30 days of HIRDLS: 3.05 GB -> 1.97 GB, for 0.65 GB of output): source
  indices are made only when `file_source=True`; no copy when the time
  selection keeps every profile; days are freed as they are stacked;
  inputs are freed once interpolated; outliers and the height range are
  masked in place. Outputs unchanged (identical in 16 real-data cases and
  the MATLAB reference tests). [Claude]
- `profiles.concat_profiles` takes `consume`, and `clip_to_height_range`
  and `mask_points` take `copy`, to allow this; defaults unchanged. [Claude]
- HIRDLS: an unknown `additional_vars` name now raises `KeyError` (it gave
  a misleading size error, or NaN with `ignore_non_modal`). [Claude]

## [0.3.3] - 2026-11-25
- `load_limb_sounders`: a clear error for `time_range` strings that are not
  real dates, e.g. "2022-01-34" (MATLAB's datenum would carry it over to
  2022-02-03): the message gives that date and `matlab_time.make_datetime`.
  [Claude]
- Test that SOFIE `Time_UT` is in hours (its files say "msec"). [Claude]
- Documentation: R3 (log-pressure interpolation) now cites eq (33b) of the
  U.S. Standard Atmosphere 1976, checked against the document, instead of
  an unchecked textbook reference; notes the small error from interpolating
  in geometric rather than geopotential height (about 5e-6 for 1 km
  levels). [Claude]
- `load_limb_sounders`, MLS: when several files match a day (e.g. v05-04
  and v06-04, as for 155 days of 2026 in the group's MLS folder), the newest
  version is now used, not the first in name order (Q8; the same change in
  the MATLAB `module_load_MLS`). Results for days with one file are
  unchanged, bit for bit.
  [Claude]
- `load_limb_sounders`: a warning, and the output attribute
  `multiple_matches`, when more than one file matches a day (any
  instrument). [Claude]
- `load_limb_sounders`: new output attributes recording the load:
  `atmosgw_curate_version`, `options` (the call's arguments, JSON),
  `files_read` and `files_not_found`; and, on each variable,
  `source_variable` (its name in the files) and a `comment` where values
  are estimated or corrected (e.g. pressure from height). Data unchanged.
  [Claude]

## [0.3.2] - 2026-11-25
- Documentation: Q5 is now fixed on the MATLAB side (`wildcardsearch.m`
  clears `regexpdir`'s persistent root before each search), so MATLAB and
  `readers.find_files` agree. No change to this package's code or results.
  [Claude]

## [0.3.1] - 2026-11-25
- `geodesy.haversine`: nearly antipodal points (within about 1e-7 degrees)
  now give half a great circle instead of NaN (Q4); the same fix in the
  MATLAB `nph_haversine` and `nph_haversine_mars`, which stopped with an error
  there. All other results unchanged, bit for bit. [Claude]

## [0.3.0] - 2026-11-25
- `atmosphere.height_to_pressure` and `pressure_to_height` now implement the
  U.S. Standard Atmosphere 1976 (Q1-Q3), as the rewritten MATLAB `h2p` and
  `p2h`: the standard's equations up to 86 km and log-pressure interpolation
  of its Table I from 86 to 149 km, valid from -5 to 149 km (NaN outside).
  Constants and table values checked against the document. The old method
  (isothermal layers between seven pressures, geopotential heights used as
  geometric) was off by up to 15% in pressure at 60 km, and 1.1 km in height.
  Effects in `load_limb_sounders`: MLS heights move (up to 0.4 km near the
  stratopause, more above), so MLS temperatures on the height grid change by
  up to 10 K; pressures computed from height change (ALOMAR, Misc and PW
  files). Other instruments are unchanged. `pressure_to_height` now always
  computes in double precision (float32 MLS pressures were computed in single
  precision). MATLAB reference outputs regenerated. [Claude]

## [0.2.0] - 2026-11-25
- `load_limb_sounders`: pressure is now interpolated to the height scale in
  log(pressure), not linearly in height (Q6, fixed the same way in the MATLAB
  `get_limbsounders`). Linear interpolation overestimated pressure between
  levels: on the MATLAB test dates by up to 11% for MIPAS (median 1.2%), 1.9%
  for MLS, 0.5% for ACE, 0.12% for HIRDLS, 0.09% for SABER and 0.02% or less
  for SOFIE, GNSS and ALOMAR. Only `pressure` changes, and only when
  interpolating (not with `original_z`, nor with `nearest`, `next` or
  `previous`). The MATLAB reference outputs were regenerated from the fixed
  MATLAB. [Claude]

## [0.1.0] - 2026-11-25
- `load_limb_sounders` checked by C. Wright, by comparing its outputs with
  the MATLAB `get_limbsounders` (SABER, HIRDLS); README and docstring
  updated to say so. [Claude]
- `load_limb_sounders(verbose=True)` now prints a line as each variable is
  interpolated, so long loads are not silent. [Claude]
- `load_limb_sounders`: the source indices are no longer interpolated when
  `file_source` is False (they were discarded afterwards), cutting
  interpolation time by about 2/7; results unchanged, bit for bit. [Claude]
- `profiles.interpolate_profiles`: linear, nearest, next and previous
  interpolation now handle all profiles at once (results unchanged, bit for
  bit); `load_limb_sounders` is now faster than the MATLAB original in every
  reference case. [Claude]
- ALOMAR: each day's own heights are now used (Q7). The first day's heights
  of a yearly file were used for every day, placing much of 2007 (from
  2007-05-14) 0.15 km too high. Fixed the same way in the MATLAB
  `module_load_ALOMAR.m`; the MATLAB reference cases are unchanged. [Claude]
- `load_limb_sounders`: added the MLS, MIPAS, SOFIE and ALOMAR readers; the
  port now covers every instrument except AIRS and IASI, and matches all 33
  MATLAB reference cases (11 bit for bit in every field; otherwise HIRDLS
  positions within 7.1e-14 degrees and pressure computed from height within
  4.9e-16 relative). [Claude]
- `atmosphere.pressure_to_height`: float32 input is computed in single
  precision, as MATLAB does (needed for MLS). [Claude]
- `load_limb_sounders`: added the HIRDLS and SABER readers. They match the
  MATLAB reference outputs for all 9 HIRDLS and SABER cases of 2007-01-01/02:
  SABER bit for bit; HIRDLS bit for bit except latitude and longitude (within
  7e-14 degrees, from the geolocation correction's trig functions). Added
  `readers.read_hdf5` and `readers.read_netcdf` (as the MATLAB `nph_getnet`).
  Adds `h5py` and `netCDF4`. [Claude]
- Added `atmosgw.curate.obs.load_limb_sounders` (port of the MATLAB
  `get_limbsounders`, in progress): the core, and the ACE/GNSS/Misc,
  planetary-wave and sampler readers. Matches the MATLAB reference outputs for
  all 14 synthetic cases (bit for bit, except pressure computed from height,
  within 4e-16 relative). The other instrument readers are still to come.
  [Claude]
- Added `atmosgw.curate.profiles` (join, select, interpolate, clip and mask
  [profiles x levels] data) and `atmosgw.curate.readers` (`find_files`,
  `read_mat`). Adds `mat-io`. [Claude]
- `matlab_compat.interp1`: linear interpolation between two equal values now
  returns exactly that value, as MATLAB does. [Claude]
- `docs/quirks.md`: added Q5 (a deliberate difference in file searching) and
  Q6 (pressure interpolated linearly in height). [Claude]
- Python 3.11 or later is now required (needed by the MATLAB file reader
  mat-io; Python 3.10 reached end of life in October 2026). The `tomli`
  dependency is no longer needed. [Claude]
- Added `atmosgw.curate.matlab_compat`: MATLAB-compatible `interp1`, `unwrap`
  and `wrap_to_pi`, checked against MATLAB R2024b output. Adds `scipy`.
  [Claude]
- Added `atmosgw.curate.atmosphere`: `height_to_pressure` and
  `pressure_to_height`, faithful ports of the MATLAB `h2p` and `p2h`, with
  their known limitations listed in `docs/quirks.md` (Q1-Q3). [Claude]
- Added `atmosgw.curate.geodesy`: spherical `haversine` (port of N. Hindley's
  `nph_haversine`), `azimuth`, `reckon` and `km2deg`, matching MATLAB.
  [Claude]
- Added `paths.dataset_dir()` and `paths.reference_dir()`: dataset folders
  default to `datasets.toml` and can be overridden in the settings file.
  [Claude]
- Added `docs/theory/relations.md` (physical relations and their sources) and
  `docs/quirks.md` (MATLAB behaviour kept for now, to fix later). [Claude]
- Added `atmosgw.curate.paths.local_data_dir()`: finds the data root from an
  argument, the `ATMOSGW_DATA` environment variable or a settings file,
  replacing hostname-based lookup. Adds `tomli` (Python 3.10 only). [Claude]
- Added `uv.lock`; CI now tests both the locked versions and the newest
  allowed versions. [Claude]
- Copyright holder in `LICENSE` is now Corwin Wright and the University of
  Bath. [Claude]
- Added `.gitattributes` pinning line endings (LF for text, CRLF for MATLAB
  and Windows scripts, no conversion for data files). [Claude]
- Added `atmosgw.curate.matlab_time`: MATLAB-compatible `datenum()` (with
  MATLAB's carry-over of out-of-range values), `make_datetime()`,
  `datenum_to_datetime()` and `datetime_to_datenum()`. [Claude]
- Added numpy and xarray (>= 2025.1.2, for non-nanosecond times) as
  dependencies. [Claude]
- Replaced AGENTS.md with CONTRIBUTING.md (contributor guidelines). [Claude]
- Local AI-agent instruction files (AGENTS.md, CLAUDE.md, .claude/) are
  git-ignored. [Claude]

## [0.0.1] - 2026-11-23
- Initial skeleton: `atmosgw` namespace package with the `atmosgw.curate`
  subpackage, packaging, licence, citation file, tests and CI. No scientific
  code. [Claude]
