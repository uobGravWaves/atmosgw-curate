# AtmosGW Curate

Python tools for loading, standardising and analysing observations of
**atmospheric gravity waves** from satellites and other instruments.

> **Status: pre-alpha (0.0.x).** Ports of the group's MATLAB tools are being
> added one at a time, and nothing should be treated as validated until it
> matches the MATLAB reference outputs. The first, `obs.load_limb_sounders`
> (from `get_limbsounders`), matches them but has not yet been reviewed.

## The AtmosGW family

AtmosGW is a small family of packages that share one import namespace:

| Package (pip name) | Import as | Purpose | Status |
|---|---|---|---|
| `atmosgw-curate` | `atmosgw.curate` | Load and standardise observations; extract gravity-wave perturbations and properties; shared data formats | This repository; skeleton |
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
    readers.py                find files; read .mat, HDF5 and netCDF files
    data/         shared data formats
    obs/          instrument loaders: load_limb_sounders (port of the MATLAB
                  get_limbsounders)
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
