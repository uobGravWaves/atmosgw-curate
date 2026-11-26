# MATLAB-vs-Python sweep of `load_limb_sounders`

A wide check of the port against the MATLAB [`get_limbsounders`](https://github.com/uobGravWaves/MultiInstrumentFunctions/blob/main/MultiInstrumentLoaders/get_limbsounders.m) on real data:
many dates for every instrument, rather than the few fixed dates of
`tests/test_limb_matlab.py`. It needs the data, MATLAB and the group's MATLAB
code, so it is run by hand (for example before a release), not by `pytest`
or CI. Nothing is written outside the work folder you give it; MATLAB outputs
are deleted as soon as they have been compared.

## Running it

```bash
export ATMOSGW_MATLAB_CODE=/path/to/MatlabFunctions  # added with genpath in MATLAB
export MATLAB=/path/to/matlab                         # default: matlab
export PYTHON=/path/to/venv/bin/python                # default: python
export SWEEP_JOBS=4                                   # MATLAB sessions at once
export SWEEP_ONLY='batch_05[0-2]'                     # optional: only these batches
tests/matlab/sweep/run_sweep.sh /path/to/empty/work/folder
```

Progress is written to `progress.log` in the work folder; a summary is printed
at the end (`summarise.py` can also be run on a work folder at any time).

## What it does

1. `make_cases.py`: for each instrument and each year with data, picks dates
   at random (fixed seed) from the days that have files: 4 per year (12 for
   HIRDLS), one of them a 3-day range, and for SABER one the first of a
   month. Each raw-data date is loaded twice: interpolated (`file_source`,
   `dx`) and with `original_z` (`file_source`). The HIRDLS, SABER, MLS and
   GNSS PW files get 2 dates per year, interpolated. Calls are grouped into
   batches of at most 30.
2. `run_batch.m`: MATLAB runs one batch of `get_limbsounders` calls and
   saves each output (or its error).
3. `compare_batch.py`: Python runs the same calls (through the private
   `_load_limb_arrays`, whose output has the MATLAB names) and compares field
   by field: field names, shapes, NaN patterns, values. One JSON line per
   call goes to `results_batch_NNN.jsonl`.
4. `summarise.py`: counts the results and explains each difference.

## Tolerances

As `tests/test_limb_matlab.py` (bit-identical unless stated):

- pressure interpolated to the height scale: 1e-14 relative (log/exp);
- pressure computed from height (ALOMAR, PW files), heights kept: 1e-15
  relative;
- HIRDLS latitude and longitude: 1e-12 degrees (azimuth/reckon trig);
- `dx`, all instruments: 1e-12 relative (haversine trig).

Differences beyond these that `summarise.py` accepts as last-bit rounding
(found in the first sweep, 2026-11-25; none affects any result):

- MLS heights computed from pressure (`p2h` / `pressure_to_height`), heights
  kept: up to 1e-15 relative (seen: 2e-16);
- pressure computed from height above 86 km, where the standard atmosphere's
  table is interpolated in log(pressure): up to 1e-14 relative (seen: 2e-15,
  ALOMAR, 97-119 km only);
- times: up to 1.2e-10 days, one unit in the last place of a datenum (10
  microseconds; outputs are rounded to the millisecond) (seen: ALOMAR, whose
  POSIX times MATLAB converts through `datetime`).

Anything else is listed as unexplained and needs looking at.

## First run

2026-11-25, MATLAB R2024b Update 6, uobGravWaves/MultiInstrumentFunctions
52518fc and corwin365/MatlabFunctions f1a78c0 (a branch made before the Q5
`wildcardsearch` fix; Q5 cannot arise here, as each MATLAB session searched
one instrument's folder with dated patterns), package 0.3.4: 1420 calls; 409 bit-identical, 807 within the tolerances above, 204
explained last-bit rounding, none unexplained; no errors; profile counts and
NaN patterns equal in every call. The six calls with no profiles are real
data gaps (SABER 2002-01-06 and 2008-01-13; two empty GNSS PW files). Total
time: MATLAB 92 min, Python 16 min.
