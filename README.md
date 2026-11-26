# AtmosGW Curate

Python tools that load observations of atmospheric gravity waves, from
satellites and ground-based instruments, into common documented formats.
Each value records why it is missing or how it was obtained, and each
Dataset records how it was loaded. Most of the loading code is ported from
the group's MATLAB code and checked against it. Wave analysis (a general
detrender, Stockwell wrappers) is planned but not yet here.

Status: pre-alpha (0.x). The data formats are version 1 proposals.

## Quick start

Requires Python 3.11 or later.

```bash
git clone https://github.com/uobGravWaves/atmosgw-curate.git
cd atmosgw-curate
uv sync --extra dev                  # exact versions from uv.lock
export ATMOSGW_DATA=/path/to/Data    # the data root: see "Where the data are"
uv run python
```

```python
from atmosgw.curate.obs import load_limb_sounders

ds = load_limb_sounders("2007-01-01", "SABER")   # one day, 18-60 km every 0.5 km
ds.temperature                                   # (profile, level), K
ds.temperature_flag                              # 0 for a plain good value; else why not
ds.attrs["files_read"]                           # the record of the load
```

`docs/usage.md` goes on from here.

## What it does today

| Function | Instruments | Returns | Origin |
|---|---|---|---|
| `obs.load_limb_sounders` | ACE, GNSS, HIRDLS, MIPAS, MLS, SABER, SOFIE, the ALOMAR lidar, Wright and Hindley sampler output, your own files in the same format; planetary-wave-filtered versions | profiles, (profile, level) | port of MATLAB `get_limbsounders` |
| `obs.load_nadir_sounders` | AIRS 3-D retrievals, AIRS 2-D brightness temperatures, IASI 3-D retrievals | swaths, (along_track, cross_track, level) | port of MATLAB `prep_airs_3d`, without its detrending and smoothing |
| `obs.load_time_height` | the ALOMAR lidar | a (time, level) series at the site | new code |
| `data.to_profiles` | swaths and time-height series | profiles | new code |

Each returns an `xarray.Dataset` in a format defined in
`docs/design/data_formats.md`. The supporting modules (MATLAB-compatible
dates and `interp1`, the U.S. Standard Atmosphere 1976, great-circle
geometry, file readers) are described in their own docstrings:
`help(atmosgw.curate.<module>)`.

## What has been checked

**`load_limb_sounders` against MATLAB `get_limbsounders`.**
`tests/test_limb_matlab.py` compares every case of a pinned MATLAB
reference file (synthetic and real data, every instrument). The test's
header names the file and the MATLAB commits; the test pins the file's
SHA-256. A by-hand sweep
at package version 0.3.4 compared 1,420 real-data calls over every
instrument and year with data (`tests/matlab/sweep/README.md`). Both compare
the private layer `_load_limb_arrays`, which has MATLAB names and datenum
times; they do not compare the public Dataset. Results are bit-identical, or
within last-bit tolerances stated in the tests. C. Wright compared the
outputs with MATLAB's on 2026-11-25 (SABER, HIRDLS); the code has not been
reviewed line by line.

**`load_nadir_sounders` against MATLAB `prep_airs_3d`.**
`tests/test_nadir_matlab.py` compares the private layer `_load_nadir_arrays`
(MATLAB-like arrays, datenum times) with MATLAB reference outputs. The test's
header names the reference and the MATLAB commit.
- AIRS 3-D, 2007-01-01: options, the dateline and midnight, on the native
  and the regridded grid; within stated last-bit tolerances.
- IASI 3-D, 2007-11-06 and 07 (two granules, load-only, midnight), with
  `screen="pooled"`, MATLAB's screen: identical wherever Python has a value.
  MATLAB refills screened points; Python leaves them NaN and flagged (under
  1 % of points in these granules). The default rolling screen differs by
  design and is not compared.
- AIRS 2-D: not compared with MATLAB. The automated tests use synthetic
  files in the files' layout; during development the core was also run on
  L. Hoffmann's real files through a throwaway reader
  (`docs/design/nadir_loader.md`, checkpoint 5).

**`load_time_height`** has no MATLAB equivalent. For ALOMAR,
`tests/test_time_height_alomar.py` compares its public Dataset with that of
`load_limb_sounders(..., "ALOMAR", original_z=True)`: the same values, rows
and states on five date ranges. A wider run on 34 dates across 2003-2026
was made once (2026-11-26).

**Metadata** (data states, `origin`, `availability`, attributes) has no
MATLAB counterpart. It is checked against the shared rules
(`docs/design/data_formats.md` section 6.3) by
`tests/test_states_conformance.py`, and by real-data audits of the state
counts. It is newer than the value checks and still being reviewed.

**Deliberate scientific changes from MATLAB**, such as log-pressure
interpolation and the 1976 standard atmosphere, are listed with their sizes
in `docs/quirks.md` and `CHANGELOG.md`. Reproducing an older result needs
the MATLAB baseline named there.

Each module's Validation section says how that module was checked.

## Where to find the detail

| Document | What it holds |
|---|---|
| `docs/usage.md` | how to load data and read what comes back |
| loader docstrings (`help(load_limb_sounders)` ...) | every option, output variable and attribute |
| `docs/design/data_formats.md` | the profile, swath and time_height formats; data states |
| `docs/design/nadir_loader.md` | design of the nadir loader and its open items |
| `docs/quirks.md` | MATLAB behaviour found wrong and fixed (Q1-Q8), with sizes |
| `docs/theory/relations.md` | physical relations [R1]-[R3], their sources and checks |
| `CHANGELOG.md` | changes by version |
| `CONTRIBUTING.md` | package layout, style, porting, documentation and tests |

## Where the data are

Loaders find instrument data under a data root folder, set once per machine
by either:

- the environment variable `ATMOSGW_DATA`, or
- a settings file `~/.config/atmosgw/config.toml` containing, for example,
  `data_dir = "C:/Data"`.

Each dataset's folder under that root is listed in
`src/atmosgw/curate/datasets.toml`. If yours are elsewhere, override just
those entries in your settings file, for example:

```toml
data_dir = "D:/Data"

[datasets.SABER]
raw = "E:/big_disk/saber"     # absolute: used as it is
pws = "SABER_v2/pws"          # relative: under data_dir
```

The full lookup rules are in `help(atmosgw.curate.paths)`.

## The AtmosGW family

AtmosGW packages share the import namespace `atmosgw`. This package,
`atmosgw-curate` (`atmosgw.curate`), loads and standardises observations
and holds everything shared: data formats, instrument characterisations and
the global-context interface. `atmosgw-fuse` (`atmosgw.fuse`, not started)
will combine curated datasets from several instruments; it depends on
Curate, never the other way round. The 2-D/3-D Stockwell transform is a
separate, general-purpose package.

## Licence, citation and contact

BSD-3-Clause; copyright Corwin Wright and the University of Bath (`LICENSE`).
To cite the package, use `CITATION.cff` (GitHub's "Cite this repository"
button). There is no DOI yet. Where a result depends on a physical
relation, also cite its source (`docs/theory/relations.md`).

Questions and bug reports: GitHub issues, or Corwin Wright, University of
Bath.

MATLAB is a registered trademark of The MathWorks, Inc. This package is
independent of MathWorks; its MATLAB-compatible functions are written from
public documentation.
