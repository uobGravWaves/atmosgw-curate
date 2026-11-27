# What has been checked

What each part of the package was compared with, at which version, and what
was not covered. Tolerances, and the largest differences seen, are in the
tests' docstrings. Runs were on the group's Linux machine (eepc-0184) with
MATLAB R2024b Update 6.

The MATLAB comparisons use the loaders' private layer (`_load_limb_arrays`,
`_load_nadir_arrays`): MATLAB names, array layout and datenum times, before
conversion to the public Dataset. The MATLAB tests (`test_get_limbsounders.m`,
`test_prep_airs_3d.m`) are in corwin365/MatlabFunctions-Tests, a private
repository.

## `load_limb_sounders`

Against MATLAB [`get_limbsounders`][get_limbsounders]:

- **`tests/test_limb_matlab.py`** (run by `pytest` when the reference files
  and data are installed; skipped otherwise, as on GitHub): 37 reference cases,
  14 synthetic and 23 real data, every instrument. Real-data dates are
  2007-01-01 and 02 for most instruments, MIPAS 2007-01-02, SOFIE
  2007-06-01, ALOMAR 2007-01-23 and 2005-01-15 (to 120 km), and MLS
  2026-01-01 (v05-04 and v06-04 files on one day). The reference was made
  2026-11-25 from MultiInstrumentFunctions@52518fc and
  MatlabFunctions@3b84a6b, which include the MATLAB fixes for Q1-Q8. Result:
  bit-identical, apart from last-bit differences in HIRDLS positions, `dx`,
  ALOMAR times, and values computed through log and exp.
- **At the port** (2026-11-24, 0.1.0), against the 33 cases then
  available (MultiInstrumentFunctions@78eef78, MatlabFunctions@c943399): 11
  cases bit-identical in every field. The rest differed only in HIRDLS
  positions (up to 7.1e-14 degrees), distances derived from them, and
  pressure computed from height (up to 4.9e-16 relative).
- **Sweep** at 0.3.4 (2026-11-25;
  [`tests/matlab/sweep/README.md`](../tests/matlab/sweep/README.md)): 1,420
  real-data calls, every instrument and year with data, interpolated and
  `original_z`, and the planetary-wave files. 409 bit-identical, 807 within
  the tolerances, 204 explained last-bit rounding, none unexplained.
- **C. Wright** compared outputs with MATLAB's (SABER, HIRDLS; 2026-11-25).
  The code has not been reviewed line by line.
- `tests/test_limb_known.py`: the MATLAB test's known-answer checks, on
  synthetic files and real data.

Each deliberate change from MATLAB (Q1-Q8) was measured before and after;
the sizes are in [`quirks.md`](quirks.md).

0.5.1 is the last release that reproduces `get_limbsounders`. Later changes
made in Python only (Q9-Q11, ACE and GNSS; Q12, boxes) are applied to the
MATLAB output in `tests/test_limb_matlab.py` before the comparison, and MLS
and SOFIE run there with `instrument_qc` off (Q15, Q17), so the rest is
still checked against MATLAB; to reproduce a MATLAB-era result, install
0.5.1.

## `load_nadir_sounders`

Against MATLAB [`prep_airs_3d`][prep_airs_3d]. IASI is run with
`screen="pooled"`, MATLAB's screen. MATLAB fills in every NaN from the
nearest point along track: screened points, and NaNs already in the files
(up to about 10 % of profiles on some days, e.g. 2007-05-31). Python leaves
them NaN and flagged, so IASI is compared wherever Python has a value.

- **`tests/test_nadir_matlab.py`** (run when the reference files and data are
  installed), reference made
  2026-11-25 from MatlabFunctions@059037e. AIRS 3-D, 8 cases on 2007-01-01:
  default, spacing, nearest, no vertical interpolation, day/night,
  load-only, the dateline and midnight. Load-only is bit-identical; the rest
  differ in the last bits. IASI 3-D, 3 cases on 2007-11-06 and 07: identical
  wherever Python has a value.
- **Sweep** at 0.4.1 (2026-11-25/26; MatlabFunctions@5112703): 878
  real-data calls. AIRS 3-D every year 2002-2025 (672 calls) and IASI 3-D
  every year 2007-2018 (206), each with default, load-only and midnight
  cases, plus dateline, day/night and spacing cases for AIRS. No
  unexplained differences. AIRS load-only bit-identical; regridded calls
  differ by at most 1.2e-10 degrees, 2.0e-12 K and 4.7e-10 days. 8 calls
  failed in both languages: 4 with no data, 4 on AIRS granules made
  entirely of fill values. The sweep's scripts and report are not in the
  repository.
- Not compared with MATLAB: the AIRS granule-file reader (every day with
  granule files also has a daily file, which both prefer), the default IASI
  screen, the `extrapolate`, `lat_range`/`lon_range`, `height_scale` and
  `n_cross`/`n_along` options, AIRS 2-D, and the fill-value handling added
  after the sweep (MATLAB reads fill values as data).

The fill-value handling (`CHANGELOG.md`, 0.5.1) comes from a census of
every AIRS daily file (8,431) and IASI file (53,439) on 2026-11-26. Against
the previous Python version, on 994 real-data calls: 974 bit-identical, 9
that crashed now load, 11 differ only at bad data.

AIRS 2-D: tests on synthetic files in the files' layout
(`tests/test_nadir_airs2d.py`); loads of real files (2020-01-20, a processed
test copy of L. Hoffmann's files) gave valid swaths. L. Hoffmann's daily
files as published (2016-01-01 and 2019-01-01, downloaded to a scratch
folder on 2026-11-27 and deleted after the checks; not kept on eepc-0184):

- **Layout:** 32,400 x 90, 240 granules of 135 rows, the first row at
  00:05:31 (granule 1); no NaN or fill values. Hour and whole-day loads
  validate; every granule's rows pass the 3-minute time check; the
  previous day's granule 240 is recorded as not found (its file was not
  downloaded); availability as for the granule files (temperature,
  pressure and winds `not_measured`).
- **Against AIRS 3-D** (the same footprints: positions equal to 4e-6
  degrees): the files' perturbations equal the brightness temperature minus
  a 4th-order polynomial fitted across track to each row, to 1e-4 K
  (4.3 and both 15 micron groups). With the same fit removed from the 3-D
  temperature at each level, the 2-D perturbation correlates best with
  the 3-D perturbation at 42 km (15 micron high, r 0.40 and 0.41) and 45 km
  (4.3 micron, r 0.26 and 0.38), with a negative lobe at 30-33 km (r down
  to -0.47), as expected of a weighting function a few scale heights deep
  over a wave of about 20 km vertical wavelength (an interpretation, not
  checked). The 15 micron low perturbation does not correlate with the 3-D
  perturbation at any height from 15 to 60 km (|r| below 0.05); not
  investigated. Peak heights from the data, not from published weighting
  functions.

CrIS 2-D (new code, 2026-11-27; L. Hoffmann's daily files):

- **Tests:** `tests/test_nadir_cris.py`, on synthetic files laid out as the
  real ones, with each field of regard's 3 x 3 pattern rotated (up to 50
  degrees) and spread towards the edges. Footprint distances from the track
  are exact to 1e-6 km; a brightness temperature linear in position is
  reproduced by the regridding to 1e-3 K; fields linear in the true along-
  and cross-track distance come back equal to the grid's own distances to
  3 km out to the swath edges, with track points zigzagging across track
  as in the real files; a missing footprint leaves at most 6 grid points
  NaN (`missing`), with every other value unchanged; the native grid is the
  index layout, footprint by footprint. Also tested: netCDF's default fill
  (9.97e36) read as missing, and granule numbers following time in a file
  with a gap.
- **Independent review** (2026-11-27, of the first version): found the
  track heading zigzagging with the 3-row scan pattern (values at the swath
  edge up to 90 km from their place along track), the fill value 9.97e36
  read as data (SNPP 2016-01-17, NOAA-20 2019-01-13), granule numbers
  counted by rows in a file with a 4.6-hour gap (SNPP 2016-01-17), and
  7.5% NaN in the outermost regridded column. All fixed; the numbers below
  are after the fixes.
- **The layout's effect, on real data** (SNPP, 2019-01-01, granules
  100-109): each regridded value lies within a median 0.2 km of its grid
  point at the centre of the swath and 1.1 km in the outer 8 columns (95%:
  4.1 and 23 km), and 99.9% of points have values. The index layout through
  the scan-line regridding used for AIRS would give 7.8 and 48.7 km (95%:
  32 and 111 km), with 86.7% of points.
- **Against the footprints** (SNPP, 2019-01-01, the whole day): regridded
  values within 4 km and 1 min of a footprint differ from it by 0.11-0.17 K
  at the 95th percentile (15 micron high, by latitude band; neighbouring
  footprints differ by up to 1.03 K); none of 322,975 by more than 2 K. The
  day's mean brightness temperature is unchanged (250.90 against 250.94 K).
  The perturbation's standard deviation falls from 0.240 to 0.190 K:
  interpolation averages footprint-to-footprint noise.
- **SNPP against NOAA-20** (2019-01-01; regridded values within 25 km and
  50 min of each other, 571,388 pairs): median difference +0.03 to +0.18 K
  (8.1 to 15 micron high), correlation 0.946 (8.1 and 10.4 micron, which
  see clouds) to 0.995 (15 micron high). 4.3 micron: 68 SNPP and 56
  NOAA-20 footprints are NaN in the files where the other groups have
  values.
- **Against AIRS 2-D** (the daily files above; AIRS points within 25 km
  and 30 min of a regridded CrIS value: about a million a day for SNPP,
  median gap 8 min; 210,000 for NOAA-20, 26 min): brightness temperatures
  correlate at 0.964-0.996, with AIRS minus CrIS medians of -1.2 K
  (4.3 micron), -1.1 K (8.1), +0.7 to +0.9 K (15 low) and -0.5 to -0.6 K
  (15 high), the channels being different averages (UNVERIFIED which).
  Perturbations correlate weakly (15 micron high 0.08-0.40, 4.3 micron
  0.14-0.63, 15 micron low near 0): the waves move in the minutes between
  the overpasses, and the footprints differ.
- **Against AIRS 3-D** (the same matching; 3-D perturbations as in the
  AIRS 2-D check): the correlation peaks at the same heights as for AIRS
  2-D (42-45 km for 15 micron high and 4.3 micron), at r 0.05-0.20.
- **Existing instruments unchanged:** 23 real-data outputs (every loader,
  `to_profiles`) bit-identical before and after, apart from the new
  `satellite` key in `options` (`null` for AIRS and IASI).
- **Not covered:** the 2022 MATLAB granules (not read); days other than
  those named; NOAA-20 on 2016-01-01 (no file).

The default IASI screen (rolling, 5 x IQR) was chosen from measurements on
four dates (`help(atmosgw.curate.obs._screen)`). A study at 0.3.6 found it
flags bad profiles (98-99 % whole profiles; T < 150 K; about 38 % exactly
100 K) and removed
no wave signal in the strongest Andes granule (2010-07-01 granule 132); a
pooled screen on a single granule removed wave crests (up to 13.5 % rms at
54 km). The study is not in the repository.

## `load_time_height`

New code, no MATLAB equivalent. `tests/test_time_height_alomar.py` compares
the public Dataset with that of `load_limb_sounders(..., "ALOMAR",
original_z=True)` on five date ranges: rows, temperature, states and
`T_ecmwf` identical; heights within 1e-9 km, times within 1 ms. The same
comparison was run once (2026-11-26) on 34 dates: a session from every
year 2003-2026, the first day of each of the 6 height grids, and an 8-day
session.

## New readers: LIMS, ROMEX, Aeolus, sampled data, trajectories

These readers are new code (2026-11-27) with no MATLAB counterpart, so they
were checked against their own files and against other data, on eepc-0184
with the package's interpreter (the scripts are not in the repository).
The ERA5 used here is 1.5 degrees and 3-hourly, so its comparisons test
units, signs, time bases and gross errors, not accuracy.

- **Tests:** `tests/test_new_loaders.py`: synthetic files with known
  answers, and the LIMS User's Guide example where the file is present.
- **Existing instruments are unchanged.** Ten limb cases (SABER, HIRDLS,
  MLS, GNSS, ACE and MIPAS on 2007-01-01; SOFIE 2007-06-01; ALOMAR
  2017-01-10; SABER with `original_z` and file sources; HIRDLS planetary
  waves) and a SWANS AIRS day are bit-identical apart from `format_version`
  (profile) and the new `options` keys. (MIPAS returned no profiles.)
- **Two days together equal the two single days**, with nothing lost or
  doubled at midnight or at the start of a file: LIMS 1978-11-08/09 and
  1979-03-26/27; ROMEX 2022-10-14/15; Aeolus 2019-01-04/05 and, for Mie,
  01-07/08; DWDsondes 2020-01-07/08 and 2004-12-31/2005-01-01; IAGOS
  2019-03-18/19; FAAM 2025-07-10/11.
- **Censuses:**
  - LIMS, all 201 files (214 days): each of the 630,053 scans is returned
    once.
  - ROMEX, every day of 2022-09 to 11: 3,126,854 profiles, every file
    profile that has a time.
  - Aeolus, the 15th of each month 2020-02 to 2023-04: every load passes
    `validate` (2022-10 and 11 are empty, a gap in the data, with a
    warning).
  - DWD: all 409 product files are in reference-time order, as the
    reader's bisection assumes.
- **Saving:** `to_netcdf`/`load_netcdf` round trips are identical unpacked,
  and within the packing resolution packed.
- **Heights**, by the hypsometric equation: LIMS `zt_std` and ROMEX
  `altitude` are geometric; LIMS `geo_heights` (to 0.04%) and DWD heights
  (to 0.06%) are geopotential.
- **ROMEX against the merged GNSS files** (2022-10-15; occultations matched
  within 1 min and about 1 degree; the 8 missions with over 100 matches):
  median difference at 18-40 km -0.26 to +0.08 K, median |difference|
  0.3-0.9 K.
- **SWANS:** the 5,552 sampled HIRDLS profiles of 2020-01-20 have the times
  of day of the real 2007-01-20 profiles, and positions within 0.06
  degrees of latitude and 0.14 of longitude (95%).
- **Against ERA5**, sampled linearly in time, bilinearly in position and
  linearly in log pressure (in height for Aeolus):
  - Aeolus HLOS (2019-01-07, 2020-07-01, 2022-06-01; both channels; 2-20
    km): with -u sin(az) - v cos(az), correlation 0.81-0.98 and slope
    0.99-1.02; with the other sign, correlation -0.81 to -0.98.
  - DWD sondes (2020-01-07, 01-08, 04-15, 07-15 and 2004-01-01):
    temperature within 0.2 K in the mean, wind slope 1.00-1.05. With the
    times shifted by 30 to 120 min either way (2020 days), the wind RMS is
    smallest, or within 0.01 m/s of the smallest, with no shift, and larger
    for every whole-hour shift. The times are therefore UTC, an offset of
    less than an hour being implausible.
  - ROMEX (2022-10-15; each mission, 10-50 km): mean difference -1.0 to
    +0.5 K, slope 0.97-1.07. COSMIC-2 heights, taken as geometric, match
    ERA5's within 30-70 m.
  - IAGOS (1998-02-27, 2010-12-04, 2019-03-18, 2023-11-15; above 400 hPa):
    0.4-0.7 K warmer than ERA5 (2023: -0.1 K), wind slope 1.01-1.03.
  - FAAM (2025-07-09, 07-15, 07-20): 0.2-0.65 K warmer; u slope 1.04-1.4
    and v correlation 0.45-0.93 (weak winds against a 1.5-degree grid).
- **Not covered:** the Aeolus profile rebuilding against the L2B product's
  own profiles (not held; the Aeolus `.mat` files do not match observation
  by observation); LIMS against independent data; other dates.
- **Code reviews.** Two independent reviews found sonde flights doubled by
  older copies in `sondes/1995/`, LIMS profiles lost where a file starts
  before its day, duplicate Aeolus observations, IAGOS flights losing their
  temperature when every package value is flagged, FAAM extras at 20 or 64
  Hz misplaced at full rate, and spurious IAGOS duplicate warnings. All
  were fixed, with tests.
- **Found in the data, and not screened:** ROMEX profiles above 330 K and
  IAGOS samples of exactly 173.16 K (`docs/data_supply.md`).
- **Later the same day** the formats (profile 3, swath 2), the choice
  between copies of a file (Q19, IAGOS, FAAM), the screening comments and
  `files_without_data` changed. Twenty-three real-data outputs (every
  loader, `to_profiles` and the new readers) were compared old against new:
  every variable is bit-identical, and only `format_version`,
  `availability` (including DWD and IAGOS `upward_wind`, now
  `not_measured`) and the new comments differ. The ERA5 comparison, rerun,
  is unchanged.

## Changes since the MATLAB comparisons

Every later change meant to leave values alone was checked old against new
on real data (3 to 39 cases each, 994 for the fill handling). Values were
identical except where intended, and times differed by up to 31
microseconds after the faster regridding (PR #25).

## Formats, states and supporting modules

Data states and the other metadata have no MATLAB counterpart.
`tests/test_data_formats.py`, `tests/test_states_conformance.py` and
`tests/test_data_convert.py` check them against the shared rules
(`docs/design/data_formats.md` section 6.3).

| Module | Checked against |
|---|---|
| `matlab_compat` | MATLAB output (`tests/data/matlab_compat_fixtures.mat`), plus 80,000 further MATLAB results for linear interpolation |
| `atmosphere` | MATLAB output of [`h2p.m`][h2p] and [`p2h.m`][p2h]; Table I of the U.S. Standard Atmosphere 1976 |
| `geodesy` | MATLAB output; `haversine` also against an independent method (10^6 pairs, within 4e-10 km) |
| `solar` | MATLAB [`which_airs_retrieval`][which_airs_retrieval] output, 2,166 places and times |
| `swath` | MATLAB output for synthetic swaths |
| `matlab_time` | MathWorks' documented examples and Python's own calendar; not MATLAB itself |
| `readers`, `paths`, `profiles` | known answers; `profiles` also through the limb MATLAB tests |
| `data.netcdf` | round trips of every format (`tests/test_netcdf.py`): unpacked identical, packed within half the stated resolution |

[get_limbsounders]: https://github.com/uobGravWaves/MultiInstrumentFunctions/blob/main/MultiInstrumentLoaders/get_limbsounders.m
[prep_airs_3d]: https://github.com/corwin365/MatlabFunctions/blob/master/DatasetSpecific/AIRS/prep_airs_3d.m
[h2p]: https://github.com/corwin365/MatlabFunctions/blob/master/GeophysicalAndTime/h2p.m
[p2h]: https://github.com/corwin365/MatlabFunctions/blob/master/GeophysicalAndTime/p2h.m
[which_airs_retrieval]: https://github.com/corwin365/MatlabFunctions/blob/master/DatasetSpecific/AIRS/which_airs_retrieval.m
