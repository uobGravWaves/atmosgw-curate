# Changelog

All notable changes to this package. Versions follow semantic versioning;
the package stays at 0.x until the Phase 1 reference tests pass.

## [Unreleased]
- Added `atmosgw.curate.matlab_time`: MATLAB-compatible `datenum()` (with
  MATLAB's carry-over of out-of-range values), `make_datetime()`,
  `datenum_to_datetime()` and `datetime_to_datenum()`. [Claude]
- Added numpy and xarray (>= 2025.1.2, for non-nanosecond times) as
  dependencies. [Claude]
- Replaced the interim AGENTS.md with the reviewed version. [Claude]

## [0.0.1] - 2026-11-23
- Initial skeleton: `atmosgw` namespace package with the `atmosgw.curate`
  subpackage, packaging, licence, citation file, tests and CI. No scientific
  code. [Claude]
