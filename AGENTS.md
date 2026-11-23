# Instructions for coding agents and contributors: atmosgw-curate

These rules apply inside this repository and take precedence over any
machine-level or personal AGENTS.md. Machine rules (environments, read-only
data directories, no `sudo`, local logging) still apply.

## What this repository is
- The Curate-stage package of the AtmosGW family (GW-LIFT / Seeing the Unseen
  programme): loaders, shared data formats, detrending, S-transform wrappers
  and wave-property calculations.
- `atmosgw` is a **PEP 420 namespace package** shared with `atmosgw-fuse` and
  later family members. **Never create `src/atmosgw/__init__.py`**; a test
  enforces this. All code for this package lives under `src/atmosgw/curate/`.
- Shared pieces (data formats, instrument characterisations, global-context
  interface) belong here. Fuse depends on Curate; this package never imports
  from Fuse.
- The S-transform itself (port of `nph_ndst`) is a **separate package and
  repository**; only gravity-wave-specific wrappers live here.
- Agree new top-level folders under `atmosgw/curate/` with the maintainers
  before adding them.

## Code style
- PEP 8 throughout, including inside functions: `snake_case` for modules,
  functions, arguments, variables and data fields; `CapWords` for classes;
  `UPPER_CASE` for constants (e.g. `N_FIXED`). This overrides any UpperCamelCase
  rule used in the group's personal MATLAB-style scripts. `ruff check .` must
  pass, including pep8-naming (`N`).
- In ports, keep the MATLAB order of operations and intermediate results,
  renamed, so lines can be compared one-to-one. Each ported function's
  docstring has a MATLAB -> Python name map.
- Functions take and return `xarray.Dataset` objects following the agreed data
  format. Times are `datetime64[us]`; use `atmosgw.curate.matlab_time` for any
  MATLAB datenum, never hand-written offsets.
- Avoid unnecessary classes and new dependencies. Adding a dependency needs the
  maintainers' approval.
- Dependencies: `pyproject.toml` states the allowed ranges; `uv.lock` records
  the exact versions used for development and CI. Change dependencies with
  `uv add` / `uv remove` (or edit `pyproject.toml` and run `uv lock`), commit
  both files together, and never edit `uv.lock` by hand. Refresh versions
  with `uv lock --upgrade` in its own commit.
- Never hard-code data paths. Use `atmosgw.curate.paths.local_data_dir()` for
  the data root; instrument sub-folders are defined by the loaders.

## MATLAB-to-Python porting traps
Check each of these in every port; all give silently different results.
- `np.interp` clamps to the end values outside the input range; MATLAB
  `interp1` returns NaN. Pass `left=np.nan, right=np.nan`.
- MATLAB is column-major: `reshape`, `(:)`, linear indices and
  `ind2sub`/`sub2ind` need `order="F"`.
- MATLAB `'` is the *conjugating* transpose; `.'` and NumPy `.T` are not.
- NumPy FFT output is not exactly conjugate-symmetric, so code that tests floats
  for exact equality (e.g. `imag == 0`) or pairs coefficients by sorting can
  break. Use index geometry instead. Never compare floats for exact equality.
- `single` vs `double`: float32 rounding can flip `max`/`argmax` ties.
- MATLAB `sort` is stable; use `kind="stable"` in NumPy.
- `wrapToPi`/`wrapTo180` differ from naive modulo at exactly +/-pi and
  +/-180 degrees.
- Indices stored in data (e.g. source file and profile numbers) stay 1-based,
  as in MATLAB, unless a decision says otherwise.

## Science rules
- Every physical relation used in code is cited by its ID from the theory and
  conventions document (`docs/theory/theory.tex`, to be added), e.g.
  `Ep from [R1]`. A relation not yet in that document is added there, marked
  UNVERIFIED, before it is used.
- Never write an equation, equation number, citation or published value from
  memory. Check it symbol by symbol against the source. If the source can't be
  accessed, mark it UNVERIFIED and say so.
- Each relation needs three checks: (P) against the source paper, (H) by an
  independent human, (N) by a numerical test. Record only checks actually done.
- Sign conventions are fixed on the conventions page of the theory document:
  upward propagation (m <= 0); cyclic wavenumbers (1/wavelength) in km^-1 in
  Phase 1; momentum flux in mPa, positive eastward/northward; direction angles
  as bearings clockwise from north, derived only. Don't introduce a new
  convention; if one is needed, stop and ask.
- Faithful port first, science fixes second. Known MATLAB bugs are reproduced
  in the port, listed in the theory document's discrepancy section, and fixed
  in a separate change with before/after tests.

## Documentation checklist (every public function or module)
1. Purpose, and where it sits in the pipeline.
2. Inputs and outputs: a table of every option (name, type, default, one-line
   description), then a breakdown of options needing more explanation.
   Shapes/dimensions, units, missing-data conventions.
3. Science: "No scientific processing" if no values change; otherwise relation
   IDs, assumptions and limits, and known caveats.
4. Origin: ported (language, repository, file, commit, authors, licence) /
   implements a published method / new method / adapted third-party code.
5. Validation: which tests, against what reference, known differences.
6. AI involvement: none / AI translation / AI-assisted / AI-originated
   (ab initio), plus who reviewed it.
7. Usage: one or two example calls.
8. Status (preliminary / operational / deprecated) and a change log of dated
   one-line entries (`YYYY/MM/DD description`); AI changes prefixed `[Claude]`.
9. How to cite: method papers, plus the package citation.

## Tests
- `ruff check .` and `pytest` must pass before work is reported as done.
  Locally: `uv sync --extra dev`, then `uv run ruff check .` and
  `uv run pytest`. CI runs both on Python 3.10 and 3.13, once with the locked
  versions and once with the newest allowed versions.
- Reference (parity) tests compare against MATLAB outputs from a pinned MATLAB
  commit. They skip when the reference files are absent and must say so; a
  skipped test is not a passed test.
- Every numerical comparison states its tolerance and why.

## Git
- Work on a feature branch and open a pull request; never push to `main`
  directly.
- Ask the maintainers before merging, tagging or releasing.
- Never commit data files, large reference outputs, credentials or virtual
  environments.
- Semantic versioning; 0.x until the Phase 1 parity tests pass. Record
  user-visible changes in `CHANGELOG.md` under "Unreleased".

## Before finishing
Summarise the files changed, tests run (interpreter, results, skips), any
relation IDs added or whose status changed, anything UNVERIFIED, and any bugs
found outside scope (reported, not fixed).
