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
| `obs.load_limb_sounders` | ACE, GNSS, HIRDLS, MIPAS, MLS, SABER, SOFIE, the ALOMAR lidar, [Wright and Hindley sampler](https://github.com/corwin365/Model-Sampling) output, your own files in the same format; planetary-wave-filtered versions | profiles, (profile, level) | port of MATLAB [`get_limbsounders`][get_limbsounders] |
| `obs.load_nadir_sounders` | AIRS 3-D retrievals, AIRS 2-D brightness temperatures, IASI 3-D retrievals | swaths, (along_track, cross_track, level) | port of MATLAB [`prep_airs_3d`][prep_airs_3d], without its detrending and smoothing |
| `obs.load_time_height` | the ALOMAR lidar | a (time, level) series at the site | new code |
| `data.to_profiles` | swaths and time-height series | profiles | new code |

Each returns an `xarray.Dataset` in a format defined in
`docs/design/data_formats.md`. The supporting modules (MATLAB-compatible
dates and `interp1`, the U.S. Standard Atmosphere 1976, great-circle
geometry, file readers) are described in their own docstrings:
`help(atmosgw.curate.<module>)`.

## What has been checked

In brief; [`docs/validation.md`](docs/validation.md) has the full record,
with versions, dates and what was not covered.

- `load_limb_sounders` matches MATLAB [`get_limbsounders`][get_limbsounders]
  in 37 reference cases (run by `pytest`) and in 1,420 real-data calls
  covering every instrument and year with data (0.3.4): bit-identical or
  last-bit differences. C. Wright has compared its outputs with MATLAB's;
  the code has not been reviewed line by line.
- `load_nadir_sounders` matches MATLAB [`prep_airs_3d`][prep_airs_3d] for
  AIRS 3-D and IASI 3-D in 11 reference cases and in 878 real-data calls
  (AIRS 2002-2025, IASI 2007-2018; 0.4.1). IASI is compared wherever Python
  has a value: MATLAB fills in missing and screened points, Python leaves
  them NaN and flagged. AIRS 2-D is not compared with MATLAB.
- These MATLAB comparisons check internal arrays (MATLAB names, datenum
  times), not the public Dataset.
- `load_time_height` (new code) gives the same values, rows and states as
  `load_limb_sounders(..., "ALOMAR", original_z=True)`.
- Data states and the other metadata have no MATLAB counterpart. They are
  checked against shared rules by tests, and are still being reviewed.
- Deliberate scientific changes from MATLAB (Q1-Q8) are listed with their
  sizes in [`docs/quirks.md`](docs/quirks.md). Reproducing an older result
  needs the MATLAB version named there.

## Where to find the detail

| Document | What it holds |
|---|---|
| `docs/usage.md` | how to load data and read what comes back |
| loader docstrings (`help(load_limb_sounders)` ...) | every option, output variable and attribute |
| `docs/design/data_formats.md` | the profile, swath and time_height formats; data states |
| `docs/validation.md` | every check against MATLAB and real data, and what was not covered |
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

AtmosGW packages share the import namespace `atmosgw`. `atmosgw-fuse`
(`atmosgw.fuse`, not started) will combine curated datasets from several
instruments. The 2-D/3-D Stockwell transform is a separate package.

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

[get_limbsounders]: https://github.com/uobGravWaves/MultiInstrumentFunctions/blob/main/MultiInstrumentLoaders/get_limbsounders.m
[prep_airs_3d]: https://github.com/corwin365/MatlabFunctions/blob/master/DatasetSpecific/AIRS/prep_airs_3d.m
