# AtmosGW Curate

Python tools that load observations of atmospheric gravity waves, from
satellites, ground-based instruments, radiosondes and aircraft, into common
documented formats.
Each value records why it is missing or how it was obtained, and each
Dataset records how it was loaded. Much of the loading code is ported from
the group's MATLAB code and checked against it; the rest is new code, checked
against other data (`docs/validation.md`). Wave analysis (a general
detrender, Stockwell wrappers) is planned but not yet here.

Status: pre-alpha (0.x). The data formats are proposals (profile version
2; swath, time_height and trajectory version 1).

## Quick start

Requires Python 3.11 or later.

```bash
git clone https://github.com/uobGravWaves/atmosgw-curate.git
cd atmosgw-curate
uv sync --extra dev                  # exact versions from uv.lock
export ATMOSGW_DATA=/path/to/Data    # the data root: see "Where the data are"
uv run python
```

Without `uv`: `pip install -e ".[dev]"` in a virtual environment. On Windows
(PowerShell), set the data root with `$env:ATMOSGW_DATA = "D:\Data"`, or
put it in `C:\Users\<you>\.config\atmosgw\config.toml`. A data folder
that does not exist is an error; days or files that are missing are listed
in `ds.attrs["files_not_found"]`, with a warning.

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
| | LIMS, ROMEX radio occultation, Aeolus HLOS winds; sampler output in netCDF (e.g. SWANS) | profiles | new code |
| `obs.load_nadir_sounders` | AIRS 3-D retrievals, AIRS 2-D brightness temperatures, IASI 3-D retrievals; model output in the AIRS layout (e.g. SWANS) | swaths, (along_track, cross_track, level) | port of MATLAB [`prep_airs_3d`][prep_airs_3d], without its detrending and smoothing |
| `obs.load_time_height` | the ALOMAR lidar | a (time, level) series at the site | new code |
| `obs.load_trajectories` | DWD high-resolution radiosondes (`DWDsondes`), IAGOS and FAAM aircraft | flights, (flight, sample) | new code |
| `data.to_profiles` | swaths and time-height series | profiles | new code |

Each returns an `xarray.Dataset` in a format defined in
`docs/design/data_formats.md`. The supporting modules (MATLAB-compatible
dates and `interp1`, the U.S. Standard Atmosphere 1976, great-circle
geometry, file readers) are described in their own docstrings:
`help(atmosgw.curate.<module>)`.

## What has been checked

0.5.1 is the last release that reproduces the MATLAB loaders. Later
releases fix problems found in review that the MATLAB loaders share (ACE and
GNSS pressure units, GNSS times and repeated occultations, box selection,
bad points blended by interpolation, nadir positions and spacing: Q9-Q14 in
[`docs/quirks.md`](docs/quirks.md)); to reproduce a MATLAB-era result,
install 0.5.1.

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
- The newer readers (LIMS, ROMEX, Aeolus, `load_trajectories`) are new code:
  known-answer tests on synthetic files, censuses of whole archives (every
  LIMS scan and ROMEX profile returned once), and comparisons with ERA5
  (temperatures, winds and the Aeolus HLOS sign) and with the merged GNSS
  files. Not yet released.
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
| `docs/data_supply.md` | the input files each loader reads, and their provenance |
| `docs/design/data_formats.md` | the profile, swath, time_height and trajectory formats; data states |
| `docs/validation.md` | every check against MATLAB and real data, and what was not covered |
| `docs/quirks.md` | MATLAB behaviour found wrong and fixed (Q1-Q18), with sizes |
| `docs/theory/relations.md` | physical relations [R1]-[R3], their sources and checks |
| `CHANGELOG.md` | changes by version |
| `CONTRIBUTING.md` | package layout, style, porting, documentation and tests |

## Where the data are

The loaders read the products below, found by file name anywhere below the
dataset's folder (`datasets.toml`): a flat folder of files works as well as
year folders. Another version of a product may work but has not been checked.
The version each file states is recorded in `ds.attrs["source_versions"]`.
`docs/data_supply.md` gives, for each product, the file names, the variables
read (with units and fill values) and where the files come from.

| Instrument | Files (as checked) | File names (usual folders) |
|---|---|---|
| ACE-FTS | v5.2 daily MAT-files | `yyyy/ace_v52_yyyydddd.mat` |
| GNSS RO | merged daily MAT-files (several centres; JPL, UCAR, ROM SAF) | `yyyy/merged_ro_yyyy_ddd.mat` |
| HIRDLS | NASA L2 v7 (HDF-EOS5) | any path containing `_yyyydddd` |
| MLS | NASA L2GP Temperature v5/v6 (HDF-EOS5); newest of several | any path containing `_yyyydddd` |
| MIPAS | Oxford retrieval v1.40 (netCDF, CEDA) | any path containing `_yyyymmdd_nom` |
| SABER | GATS custom Level 2A monthly Temp_O3_H2O files (v2.07, v2.08) | any path containing `<Month><yyyy>` |
| SOFIE | GATS Level 2 v1.3 (netCDF) | any path containing `Level2_yyyyddd` |
| ALOMAR lidar | yearly HDF5, V2026-07-27 | `yyyy_V2026-07-27_file4juli.h5` |
| AIRS 3-D | L. Hoffmann's retrieval: daily files, else granule files | `yyyy/airs_yyyydddd.nc`, `yyyy/ddd/airs_yyyy_ddd_ggg.nc` |
| AIRS 2-D | L. Hoffmann's brightness-temperature granule files | `yyyy/ddd/airs_yyyy_ddd_ggg.nc` |
| IASI 3-D | L. Hoffmann's retrieval, orbit files | `yyyy/ddd/iasi_yyyymmddhhmmss.nc` |
| LIMS | V6 revised Level 2 ASCII | `LIMS_2_yyyymmdd_V6.ASCII` |
| ROMEX | GNSS RO per-mission netCDF, v1a | `romex_<mission>rmx..._100m_v1a....nc` |
| Aeolus | L2B HLOS winds, monthly flat netCDF (VirES) | `Aeolus_L2B_monthly_yyyy_mm.nc` |
| DWD radiosondes | high-resolution ("Sekundenwerte aero") text files | `produkt_sec_aero_yyyymmdd_yyyymmdd_<station>.txt` |
| IAGOS | L2 time series, one netCDF per flight (in a zip, or unzipped) | `IAGOS_timeseries_<yyyymmddHHMMSS><nn>....nc[4]` |
| FAAM | core data, one netCDF per flight (v005) | `core_faam_yyyymmdd_v005_r<rev>_<flight>.nc` |
| SWANS (model) | sampler netCDF output; AIRS granule layout | `<INST>-dyamond_<model>-yyyydddd.nc`; as AIRS 3-D |

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
Bath. With a report, please paste the output's `ds.attrs` (the call's
options, the package version, the files read and their data versions):
with them, most problems can be reproduced.

MATLAB is a registered trademark of The MathWorks, Inc. This package is
independent of MathWorks; its MATLAB-compatible functions are written from
public documentation.

[get_limbsounders]: https://github.com/uobGravWaves/MultiInstrumentFunctions/blob/main/MultiInstrumentLoaders/get_limbsounders.m
[prep_airs_3d]: https://github.com/corwin365/MatlabFunctions/blob/master/DatasetSpecific/AIRS/prep_airs_3d.m
