# Contributing to AtmosGW Curate

## Package layout
- `atmosgw` is a namespace package shared with `atmosgw-fuse` and later family
  members. Never add `src/atmosgw/__init__.py`; a test checks this. All code
  for this package lives under `src/atmosgw/curate/`.
- Shared pieces (data formats, instrument characterisations, global-context
  interface) belong here. `atmosgw-fuse` depends on this package, not the
  other way round.
- The Stockwell transform is a separate package; only gravity-wave-specific
  wrappers live here.
- Agree new top-level folders under `atmosgw/curate/` with the maintainers
  first.

## Style
- PEP 8 naming throughout, including inside functions (`snake_case`
  variables and functions, `CapWords` classes, `UPPER_CASE` constants).
  `ruff check .` enforces this.
- Functions take and return `xarray.Dataset` objects. Times are
  `datetime64[us]`; convert MATLAB datenums with `atmosgw.curate.matlab_time`.
- Find data with `atmosgw.curate.paths.local_data_dir()`; never hard-code
  paths.
- Keep dependencies few. Change them with `uv add` / `uv remove` and commit
  `pyproject.toml` and `uv.lock` together.

## Porting from MATLAB
- Match the MATLAB reference outputs (from a pinned MATLAB commit) before
  changing any behaviour. Science fixes come afterwards, each in its own
  change with before/after tests.
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
- Physical relations are listed, with sources, in the theory and conventions
  document (`docs/theory/`), and cited in code by ID, e.g. `[R2]`.
- Equations and citations are checked against the original paper, by a
  second person, and by a numerical test. Anything not yet checked is marked
  UNVERIFIED.
- Conventions: upward propagation (m ≤ 0); wavenumbers as 1/wavelength in
  km⁻¹; momentum flux in mPa, positive eastward/northward; directions as
  bearings clockwise from north.

## Documentation
Each public function or module states: purpose; inputs and outputs (a table
of options, then details); science (or "No scientific processing"); origin;
validation; AI involvement; a usage example; status and a dated change log;
how to cite.

## Tests
- `uv sync --extra dev`, then `uv run ruff check .` and `uv run pytest`.
  CI runs both on Python 3.10 and 3.13, with locked and with latest versions.
- MATLAB reference tests skip when their data are absent; a skip is not a
  pass.
- State the tolerance in every numerical comparison.

## Git
- Work on a branch and open a pull request; don't push to `main`.
- Don't commit data, large reference outputs or environments.
- Semantic versioning (0.x until the ports match MATLAB). Note
  user-visible changes in `CHANGELOG.md` under "Unreleased".
