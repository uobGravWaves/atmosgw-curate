# AtmosGW Curate

Python tools for loading, standardising and analysing observations of
**atmospheric gravity waves** from satellites and other instruments.

> **Status: pre-alpha (0.x), with three working loaders.**
> - `obs.load_limb_sounders` (port of `get_limbsounders`): ACE, GNSS, HIRDLS,
>   MIPAS, MLS, SABER, SOFIE, the ALOMAR lidar, sampler output and files of
>   your own, as profiles.
> - `obs.load_nadir_sounders` (port of `prep_airs_3d`, without its
>   detrending): AIRS 3-D retrievals, AIRS 2-D brightness temperatures and
>   IASI 3-D retrievals, as swaths (`data.to_profiles` turns them into
>   profiles).
> - `obs.load_time_height` (new code): the ALOMAR lidar as a [time x level]
>   series at its site (also convertible to profiles).
>
> All return the common formats of `docs/design/data_formats.md` (version 1,
> a proposal). See `docs/usage.md` to get started.

## Validation status

- **Values against MATLAB.** `load_limb_sounders` matches the MATLAB
  `get_limbsounders` in its regression cases and in a 1,420-call real-data
  sweep at 0.3.4 (`tests/matlab/sweep/README.md`); checked by C. Wright on
  2026-11-25 (SABER, HIRDLS). `load_nadir_sounders` matches `prep_airs_3d` for
  AIRS (dateline, midnight, native and regridded) and, for IASI, wherever
  Python keeps a value (its default screen differs by design).
  `load_time_height` has no MATLAB equivalent: for ALOMAR it gives the same
  values and rows as `load_limb_sounders(..., original_z=True)` on 34 dates
  across 2003-2026 (`tests/test_time_height_alomar.py`). Each reference
  file is pinned to the MATLAB commits named in its test's header
  (`tests/test_limb_matlab.py`, `tests/test_nadir_matlab.py`).
- **Metadata (data states, `origin`, `availability`, attributes).** MATLAB has
  none of this, so it is checked differently: by the shared rules in
  `docs/design/data_formats.md` section 6.3, the cross-loader tests
  `tests/test_states_conformance.py`, and real-data audits of the state counts.
  It is newer than the value checks and still being reviewed.
- **Scientific changes from MATLAB** (e.g. log-pressure interpolation, the
  1976 standard atmosphere) are listed with their sizes in `CHANGELOG.md` and
  `docs/quirks.md`; reproducing an older result needs the MATLAB baseline
  named there.

## The AtmosGW family

AtmosGW is a small family of packages that share one import namespace:

| Package (pip name) | Import as | Purpose | Status |
|---|---|---|---|
| `atmosgw-curate` | `atmosgw.curate` | Load and standardise observations; extract gravity-wave perturbations and properties; shared data formats | This repository; loaders working, analysis to come |
| `atmosgw-fuse` | `atmosgw.fuse` | Combine curated datasets from multiple instruments | Not started |

Curate holds everything shared (data formats, instrument characterisations,
global-context interface). Fuse depends on Curate, never the other way round.

The 2-D/3-D Stockwell transform used for wave analysis is a separate,
general-purpose package outside this family.

## Installation (development)

```bash
git clone https://github.com/uobGravWaves/atmosgw-curate.git
cd atmosgw-curate
uv sync --extra dev        # exact versions from uv.lock
uv run pytest
uv run ruff check .
```

`pip install -e ".[dev]"` also works, with the newest allowed versions
instead of the locked ones. Requires Python 3.11 or later.

## Using it

See `docs/usage.md` for loading limb profiles (`obs.load_limb_sounders`) and
nadir swaths (`obs.load_nadir_sounders`), and reading what they return.

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

Tests that compare against MATLAB reference outputs look for them in
`reference_dir` (set in the same file; relative paths are under `data_dir`),
and skip if it is not set.

## Planned layout

Agree changes with the maintainers before adding to this. Folders are added
only when code for them exists.

```
src/atmosgw/curate/
    paths.py, datasets.toml   where the data are
    matlab_time.py            MATLAB-compatible dates
    matlab_compat.py          MATLAB-compatible interp1, unwrap, wrapToPi
    atmosphere.py             height <-> pressure (approximate)
    geodesy.py                great-circle distance, azimuth, reckon
    profiles.py               join, select, interpolate, clip and mask profiles
    swath.py                  distances, regridding and vertical grids of swaths
    solar.py                  solar zenith angle; AIRS/IASI day-night switch
    readers.py                find files; read .mat, HDF5 and netCDF files
    data/         shared data formats: swath and profile formats, data
                  states, validate, to_profiles (proposal)
    obs/          instrument loaders: load_limb_sounders (port of the MATLAB
                  get_limbsounders); load_nadir_sounders (AIRS, IASI; port
                  of prep_airs_3d; preliminary)
    instrument/   instrument characterisations (CPs)
    context/      global context (GC) interface
    analysis/     detrending, Stockwell wrappers
    retrieval/    parameter engine
    sampling/     sampler and synthetic observations
```

**Important:** `atmosgw` is a namespace package. Never add
`src/atmosgw/__init__.py`; doing so would stop `atmosgw-curate` and
`atmosgw-fuse` being installed together. A test enforces this.

## Conventions

- PEP 8 naming throughout, including inside functions; enforced by `ruff`
  (rule set `N`).
- Every function follows the documentation checklist in `CONTRIBUTING.md`.
- Semantic versioning; stays at 0.x until the ports match the MATLAB
  reference outputs.

## Licence and citation

BSD-3-Clause; copyright Corwin Wright and the University of Bath. See `LICENSE`. Citation details are in `CITATION.cff`
(GitHub's "Cite this repository" button). A Zenodo DOI will be added with the
first release.

## Trademarks

MATLAB is a registered trademark of The MathWorks, Inc. This package is
independent of MathWorks; its MATLAB-compatible functions are written from
public documentation.

## Contact

Corwin Wright, University of Bath. Please use GitHub issues for questions and
bug reports.
