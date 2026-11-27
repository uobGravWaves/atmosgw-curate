# Contributing to AtmosGW Curate

## Package layout
- `atmosgw` is a namespace package shared with `atmosgw-fuse` and later family
  members. Never add `src/atmosgw/__init__.py`: it would stop the packages
  being installed together, and a test checks for it. All code for this
  package lives under `src/atmosgw/curate/`.
- Shared pieces (data formats, instrument characterisations, global-context
  interface) belong here. `atmosgw-fuse` depends on this package, not the
  other way round.
- The Stockwell transform is a separate package; only gravity-wave-specific
  wrappers live here.
- Agree new top-level folders under `atmosgw/curate/` with the maintainers
  first. Folders are added only when code for them exists. The planned
  layout:

```
src/atmosgw/curate/
    paths.py, datasets.toml   where the data are
    matlab_time.py            MATLAB-compatible dates
    matlab_compat.py          MATLAB-compatible interp1, unwrap, wrapToPi
    atmosphere.py             height <-> pressure (U.S. Standard Atmosphere 1976)
    geodesy.py                great-circle distance, azimuth, reckon
    profiles.py               join, select, interpolate, clip and mask profiles
    swath.py                  distances, regridding and vertical grids of swaths
    solar.py                  solar zenith angle; AIRS/IASI day-night switch
    readers.py                find files; read .mat, HDF5 and netCDF files
    data/         data formats, data states, validate, to_profiles (proposal)
    obs/          instrument loaders: load_limb_sounders, load_nadir_sounders
                  (preliminary), load_time_height, load_trajectories
    instrument/   instrument characterisations (CPs)        (planned)
    context/      global context (GC) interface             (planned)
    analysis/     detrending, Stockwell wrappers            (planned)
    retrieval/    parameter engine                          (planned)
    sampling/     sampler and synthetic observations        (planned)
```

## Style
- PEP 8 naming throughout, including inside functions (`snake_case`
  variables and functions, `CapWords` classes, `UPPER_CASE` constants).
  `ruff check .` enforces this (rule set `N`).
- Functions take and return `xarray.Dataset` objects. Times are
  `datetime64[us]`; convert MATLAB datenums with `atmosgw.curate.matlab_time`.
- Find data with `atmosgw.curate.paths.local_data_dir()`; never hard-code
  paths.
- Keep dependencies few. Change them with `uv add` / `uv remove` and commit
  `pyproject.toml` and `uv.lock` together.

## Porting from MATLAB
- Match the MATLAB reference outputs (from a pinned MATLAB commit) before
  changing any behaviour. Science fixes come afterwards, each in its own
  change with before/after tests, and are listed in `docs/quirks.md`.
- Keep the MATLAB order of operations so code can be compared line by line,
  and give a MATLAB → Python name map in the docstring.
- Common traps: `np.interp` clamps where `interp1` returns NaN (pass
  `left=np.nan, right=np.nan`); MATLAB is column-major (`order="F"`); MATLAB
  `'` is a conjugating transpose; `sort` is stable in MATLAB
  (`kind="stable"`); float32 can flip `argmax` ties; never test floats for
  exact equality; stored source indices stay 1-based.

## Code from other sources
- **MathWorks functions** (anything under `matlabroot`, e.g. `wrapToPi`,
  `interp1`): never copy or translate their source or documentation.
  Re-implement from the documented behaviour and test against MATLAB outputs.
- **Third-party code** (e.g. File Exchange): port only if the licence allows,
  and keep its notice.
- **The group's own code:** port freely.

Record which applies in the docstring's Origin section.

## Science
- Physical relations are listed, with sources, in
  `docs/theory/relations.md`, and cited in code by ID, e.g. `[R2]`.
- Equations and citations are checked against the original paper, by a
  second person, and by a numerical test. Anything not yet checked is marked
  UNVERIFIED.
- Conventions: upward propagation (m ≤ 0); wavenumbers as 1/wavelength in
  km⁻¹; momentum flux in mPa, positive eastward/northward; directions as
  bearings clockwise from north.

## Documentation
Each fact has one home; elsewhere, link to it. `README.md` says what the
package does and summarises what has been checked; `docs/validation.md` is
the full record of checks (add each new comparison with MATLAB or real
data there); `docs/usage.md` shows common tasks; the module docstrings are
the reference for options, outputs and science.

Module docstrings are read through `help()`, often without the repository
open, so each should make sense on its own. Start with a one-line summary
and a paragraph on what the module is for. Then use these numpydoc
sections, leaving out any with nothing specific to say:

- **Inputs and outputs**: options (a table, or a list with type and default),
  then details: units, dimensions and their order, defaults, missing-data
  behaviour, and any MATLAB behaviour reproduced or deliberately changed.
- **Science**: equations and assumptions, each with its source or relation
  ID. Mark anything unchecked UNVERIFIED.
- **Origin**: the MATLAB file and commit ported, or "new code", and which
  rule under "Code from other sources" applies.
- **Validation**: what was compared with what (reference version, internal
  arrays or the public output), and where the tolerances are stated.
- **Status**: preliminary, proposal or matching MATLAB; who wrote it (AI
  involvement, e.g. "AI-ported (Claude)") and who has reviewed it.
- **Usage**: a short example.
- **Change log**: one dated line per change, e.g.
  `2026/11/25 [Claude] pressure interpolated in log(pressure) (Q6)`.

Private modules and simple helpers can be much shorter: a few sentences on
what they do, the file layout they read, and their differences from MATLAB.
Public functions whose module documents them can say so in one line.
Citation instructions belong in `README.md` and `CITATION.cff`, not in each
module.

## Tests
- `uv sync --extra dev`, then `uv run ruff check .` and `uv run pytest`.
  `pip install -e ".[dev]"` also works, with the newest allowed versions
  instead of the locked ones. CI runs ruff and pytest on Python 3.11 and 3.13,
with locked and with latest versions.
- Tests that compare against MATLAB reference outputs find them in
  `reference_dir`, set in the settings file of `README.md` ("Where the data
  are"; relative paths are under `data_dir`). They skip if it is not set or
  the data are absent; a skip is not a pass.
- State the tolerance in every numerical comparison.

## Git
- Work on a branch and open a pull request; don't push to `main`.
- Don't commit data, large reference outputs or environments.
- Semantic versioning (0.x until the ports match MATLAB). Note
  user-visible changes in `CHANGELOG.md` under "Unreleased".
