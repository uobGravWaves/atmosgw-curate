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
made in Python only (Q9-Q11, ACE and GNSS) are applied to the MATLAB output
in `tests/test_limb_matlab.py` before the comparison, so the rest is still
checked against MATLAB; to reproduce a MATLAB-era result, install 0.5.1.

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

The fill-value handling (`CHANGELOG.md`, Unreleased) comes from a census of
every AIRS daily file (8,431) and IASI file (53,439) on 2026-11-26. Against
the previous Python version, on 994 real-data calls: 974 bit-identical, 9
that crashed now load, 11 differ only at bad data.

AIRS 2-D: tests on synthetic files in the files' layout
(`tests/test_nadir_airs2d.py`); loads of real files (2020-01-20, a processed
test copy of L. Hoffmann's files; the originals are not on eepc-0184) gave
valid swaths.

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

New code (2026-11-27), with no MATLAB counterpart: LIMS, ROMEX and Aeolus
(`load_limb_sounders`), netCDF sampler output and `sampled_from` (limb and
nadir), and `load_trajectories` (DWDsondes, IAGOS, FAAM). Checked as below,
with the package's own interpreter on eepc-0184; the scripts are not in the
repository. Differences from ERA5 include ERA5's coarse resolution (1.5
degrees, 3-hourly, as kept on eepc-0184) and real waves, so they test units,
signs, time bases and gross errors, not the data's accuracy.

- **Tests:** `tests/test_new_loaders.py`, synthetic files of each format
  (known answers for parsing, units, times, fill values, flags, profile
  rebuilding, file choice, options), plus the User's Guide example of a
  real LIMS file where present.
- **Existing instruments unchanged:** ten limb cases (SABER, HIRDLS, MLS,
  GNSS, ACE, MIPAS 2007-01-01; SOFIE 2007-06-01; ALOMAR 2017-01-10; SABER
  with original heights and file sources; HIRDLS planetary-wave data) and a
  SWANS AIRS day, before and after: every variable and variable attribute
  bit-identical; only `format_version` (profile 1 to 2) and the new options
  in `options` differ.
- **One request over two days equals the two single days** (no profile or
  flight lost or doubled at midnight or at a file's start): LIMS
  1978-11-08/09 and 1979-03-26/27, ROMEX 2022-10-14/15, Aeolus 2019-01-04/05
  (across a file start) and 2019-01-07/08 (Mie), DWDsondes 2020-01-07/08
  and 2004-12-31/2005-01-01, IAGOS 2019-03-18/19, FAAM 2025-07-10/11.
- **Censuses:** LIMS, every day of the 201 files: the 630,053 scans are
  each returned once, and every day validates. ROMEX, every day
  2022-09-01 to 11-30: 3,126,854 profiles returned, exactly the files'
  profiles with a time. Aeolus, the 15th of every month 2020-02 to 2023-04:
  every day validates; the files have no data from 2022-10-05 to 2022-11-26
  (a gap in the files: those loads are empty, with a warning). DWD sondes:
  all 409 product files are in order of reference time, as the reader's
  bisection assumes (it also falls back to a full read if not).
- **Saving:** `data.to_netcdf`/`load_netcdf` round trips of each new
  instrument's output: identical, and within 0.01 when packed.
- **Heights against hypsometric thickness:** LIMS `zt_std` is geometric
  (increments exceed the thickness by the geometric/geopotential ratio,
  0.5-1.5% at 10-50 km) and `geo_heights` geopotential (to 0.04%); ROMEX
  altitude is geometric (as LIMS); DWD heights are geopotential (to 0.06%).
- **ROMEX against the merged GNSS files** (2022-10-15, the same occultations
  matched within 1 min and about 1 degree): for the eight missions with
  over 100 matches (COSMIC-2, KOMPSAT-5, MetOp-B, MetOp-C, PAZ, Spire,
  TanDEM-X, TerraSAR-X), median difference at 18-40 km -0.26 to +0.08 K,
  median |difference| 0.3-0.9 K. FengYun-3, Sentinel-6, Tianmu, Yunyao and
  GeoOptics have few or no counterparts in the merged files.
- **SWANS:** the sampled HIRDLS profiles of 2020-01-20 have the times of day
  of the real HIRDLS profiles of 2007-01-20 they were sampled from (all
  5,552 matched by profile number), and positions within 0.06 degrees
  latitude and 0.14 degrees longitude (95%; mean over levels).
- **Against ERA5** (local daily model-level files; linear in time, bilinear
  in position, linear in log pressure, or in geopotential height for
  Aeolus):
  - Aeolus HLOS, 2019-01-07, 2020-07-01, 2022-06-01, both channels, 2-20
    km: with HLOS = -u sin(az) - v cos(az), correlation 0.81-0.98 and slope
    0.99-1.02; with the other sign -0.81 to -0.98. This is the sign
    convention now stated in `hlos_wind`'s `sign_convention`.
  - DWD sondes, four days of 2020 and 2004-01-01: temperature within 0.2 K
    (mean), winds slope 1.00-1.05; and the wind differences are smallest
    with the times as read, growing for shifts of 30, 60 and 120 min either
    way, so the files' times are UTC.
  - ROMEX, 2022-10-15, every mission 10-50 km: mean differences -1.0 to
    +0.5 K, slope 0.97-1.07. COSMIC-2 heights, taken as geometric, match
    ERA5's at the same pressure within 30-70 m.
  - IAGOS, flights of 1998-02-27, 2010-12-04, 2019-03-18 and 2023-11-15
    above 400 hPa: temperature 0.4-0.7 K warmer than ERA5 (except 2023:
    -0.1 K), winds slope 1.01-1.03.
  - FAAM, 2025-07-09, 07-15, 07-20: temperature 0.2-0.65 K warmer, u slope
    1.04-1.4 and v correlation 0.45-0.93 (weak winds at the flights'
    levels, against ERA5's 1.5-degree grid).
- **Not covered:** the Aeolus profile rebuilding against the L2B
  product's own profiles (not held on eepc-0184); the separate Aeolus `.mat`
  extraction (2018-2022) holds other observations or processing (no
  observation-by-observation agreement), so it could not serve as a check;
  LIMS against independent data (no reanalysis for 1978-79 on eepc-0184);
  DWD sondes beyond the dates above (1993-06-15, 2004-01-01, 2020, 2024-12-31
  were loaded); IAGOS beyond four days (1998, 2010, 2019, 2023) of its 30
  years; FAAM beyond its 27 flights of 2025-07-09 to 08-01.
- **Independent reviews** (one subagent each, of the code without the
  author's reasoning): the first found sonde flights doubled by the older
  file copies in `sondes/1995/`, LIMS profiles lost where a file starts
  before its day, and duplicate Aeolus observations; the second found IAGOS
  flights losing their temperature when the package values are all
  flagged, FAAM extras at 20 or 64 Hz misplaced at full rate, and spurious
  IAGOS duplicate warnings. All were fixed, with tests, and the real-data
  cases rechecked.
- **Found in the data (reported, not screened):** 22-45 ROMEX profiles a day
  (about 0.1%) not marked `is_bad` exceed 330 K above 30 km; one 2023 IAGOS
  flight has 82 take-off samples of exactly 173.16 K marked good.

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
