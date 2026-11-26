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

[get_limbsounders]: https://github.com/uobGravWaves/MultiInstrumentFunctions/blob/main/MultiInstrumentLoaders/get_limbsounders.m
[prep_airs_3d]: https://github.com/corwin365/MatlabFunctions/blob/master/DatasetSpecific/AIRS/prep_airs_3d.m
[h2p]: https://github.com/corwin365/MatlabFunctions/blob/master/GeophysicalAndTime/h2p.m
[p2h]: https://github.com/corwin365/MatlabFunctions/blob/master/GeophysicalAndTime/p2h.m
[which_airs_retrieval]: https://github.com/corwin365/MatlabFunctions/blob/master/DatasetSpecific/AIRS/which_airs_retrieval.m
