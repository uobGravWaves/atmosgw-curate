# Notes for coding agents and contributors

## Package structure
- `atmosgw` is a **PEP 420 namespace package** shared by `atmosgw-curate` (this
  repository) and `atmosgw-fuse`. **Never create `src/atmosgw/__init__.py`.**
  All code for this package lives under `src/atmosgw/curate/`.
- Shared pieces (data formats, instrument characterisations, global-context
  interface) belong here in Curate. Fuse depends on Curate, never the reverse.
- Agree new top-level folders under `atmosgw/curate/` with the maintainers
  before adding them.

## Style
- PEP 8 naming everywhere, including inside function bodies: `snake_case` for
  modules, functions, arguments, variables and data fields; `CapWords` for
  classes; `UPPER_CASE` for constants. This overrides any UpperCamelCase rule
  used in the group's personal MATLAB-style scripts.
- Run `ruff check .` and `pytest` before committing.

## Porting from MATLAB
- Faithful port first, science fixes second: match the MATLAB reference
  outputs (pinned to a specific MATLAB commit) before changing behaviour.
- Keep a MATLAB → Python name map in each ported function's docstring.
- Known traps: `np.interp` clamps where MATLAB returns NaN; MATLAB is
  column-major (use `order="F"`); MATLAB `'` is a conjugating transpose; never
  compare floats for exact equality.

## Documentation
Every function follows the group documentation checklist: purpose; inputs and
outputs; science; origin; validation; AI involvement; usage; status and
changes; how to cite. Tag AI-made changes in the function's change log
(e.g. `[Claude]`).
