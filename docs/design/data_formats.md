# Design: common data formats

**Status:** draft for agreement, 2026-11-25. Nothing here is implemented yet.
Points marked *agreed* were agreed with C. Wright on 2026-11-25.

The Curate plan gives "common interfaces and data structures" very high
priority: every loader's output must be usable by the same Detrender,
Stockwell wrappers, Parameter Engine, sampler and Fuse. This document defines
those outputs as named, versioned **formats**, so that "parallel interfaces"
can be checked, not just agreed.

Formats now: **profile** (limb sounders, occultation, lidar; the output of
`load_limb_sounders`) and **swath** (nadir sounders; `load_nadir_sounders`,
see `nadir_loader.md`). Later, as loaders for them are written:
**time-height** (radars, lidar time series), **point series** (balloons) and
**image** (ground-based imagers).

The specification will live in code in `atmosgw/curate/data/` (a new folder,
already in the README's planned layout): `formats.py` (the definitions below
as data), `states.py` (the registry of data states, §6) and
`validate(ds, format)`.

## 1. What every format specifies

1. **Dimensions**: names, order, which are optional, and the length-0 case.
2. **Required variables**: name, dimensions, type, units, CF `standard_name`
   or `long_name`.
3. **Optional variables**: by name or naming pattern.
4. **Required dataset attributes.**
5. **Dimensionality metadata**: which dimensions are spatial or temporal, and
   their spacing, for the Detrender and the CP's dimensionality constraints.
6. **Conversions** to other formats.
7. **Version**: an integer, raised on any incompatible change.

## 2. Conventions common to all formats

- Times: `datetime64[us]`, UTC, named `time`. Missing: NaT. Leap seconds are
  ignored (*agreed*: well below wave timescales).
- Latitude: degrees north, [-90, 90]. Longitude: degrees east, [-180, 180].
- Altitude: km (geometric, as given by each instrument's data; a variable
  `comment` says if it is geopotential). Pressure: hPa. Temperature: K.
- No value: NaN (NaT for times), always with a state saying why (§6). Bad
  points are NaN in every loader, never refilled (*agreed*). No fill values
  in outputs.
- Source indices (rows, profiles, files): 1-based, as stored in the data.
- Signs and wave conventions: as `CONTRIBUTING.md` (upward propagation
  m <= 0; wavenumbers as 1/wavelength in km^-1; bearings clockwise from
  north).
- Names: lower case, words separated by `_`, spelled out (`temperature`, not
  `T`).
- Attributes are text, numbers or lists of them, so a Dataset survives
  `to_netcdf` unchanged.

### Numeric types

- **Data variables: float32.** Its resolution is about 1 part in 1.7e7:
  1.5e-5 K at 300 K, 1e-5 degrees (about 1 m) in longitude. That is far finer
  than the data (about 1 K for nadir sounders, about 0.1 K at best for limb
  sounders and GNSS), so nothing is lost, and memory is half of float64.
  *Proposal (2026-11-25):* the "profile" format also accepts float64 data
  variables, so `load_limb_sounders` keeps its values exactly as loaded
  (most of its readers give float64); the validator checks this per format.
- **Coordinates and distances: float64**; times `datetime64[us]`.
- **Calculations: float64.** Analysis functions convert on input and return
  float32, so rounding never accumulates in float32.
- **State flags:** unsigned integers (§6).
- **Packing (scale and offset):** only for files written to disk, never in
  memory (xarray unpacks on reading). See §11.

## 3. Format "profile", version 1

**Dimensions:** `profile`, `level`; optional `mode` (planetary-wave modes),
`file`. `profile` may have length 0.

| Variable | Dimensions | Units | Required |
|---|---|---|---|
| `latitude`, `longitude` | (profile, level) | degrees | yes |
| `time` | (profile, level) | — | yes |
| `altitude` | (profile, level) | km | yes |
| `temperature` | (profile, level) | K | yes, unless PW data |
| `pressure` | (profile, level) | hPa | yes; `origin` attribute says whether measured or derived (§6.3) |
| `temperature_residual`, `temperature_pw` | (profile, level[, mode]) | K | PW data |
| `<variable>_flag` | as `<variable>` | — | yes, for each data variable (§6) |
| `source_profile`, `source_file` | (profile, level) | — | optional |
| `original_files` | (file) | — | optional |
| `dx` | (profile, level) | km | optional |

This is the `load_limb_sounders` output since 2026-11-25 (values unchanged by
the format; data variables may be float32 or float64). States it sets:
`bad_geolocation`, `outside_physical_range`, `outside_requested_range` (with
`original_z`),
`no_observation` (every other NaN), `extrapolated` (with `extrapolate`),
`derived` (estimated pressure, point by point;
planetary-wave data, Sampler data). *Proposals:* `product` "temperature" and
`characterisation` `<instrument>-TEMPERATURE-placeholder/0` for every limb
instrument; `origin` "retrieved" for limb temperatures and pressures.

## 4. Format "swath", version 1

**Dimensions:** `along_track`, `cross_track`, then `level` (3-D products) or
`channel` (2-D products with several channels; absent with one). The order is
fixed. `along_track` may have length 0.

| Variable | Dimensions | Units | Required |
|---|---|---|---|
| `latitude`, `longitude` | (along_track, cross_track) | degrees | yes |
| `time` | (along_track, cross_track) | — | yes |
| `altitude` | (level) | km | 3-D; on the native grid, (along_track, cross_track, level) if heights vary by point |
| `temperature` | (along_track, cross_track, level) | K | 3-D |
| `brightness_temperature` | (along_track, cross_track[, channel]) | K | 2-D |
| `channel_wavenumber` | (channel) | cm^-1 | 2-D with channels |
| `channel_name` | (channel) | — | 2-D with channels: the channel or channel group |
| `brightness_temperature_perturbation`, `brightness_temperature_variance` | as `brightness_temperature` | K, K2 | optional: values derived by the files' producer (`origin` "derived") |
| `<variable>_flag` | as `<variable>` | — | yes, for each data variable (§6) |
| `granule`, `granule_day` | (along_track) | — | yes: the granule each row came from, and its day |
| `segment` | (along_track) | — | yes: 1-based number of the unbroken stretch the row is in |
| `distance_along`, `distance_cross` | (along_track), (cross_track) | km | yes on a regular grid |
| `source_file`, `source_row` | (along_track) | — | optional |
| `original_files` | (file) | — | optional |
| `day_retrieval` | (along_track, cross_track) | — | optional |
| native geolocation: `native_latitude`, `native_longitude`, `native_time`, `scan_position` | (native_along, native_cross) | — | optional |
| viewing geometry: `satellite_latitude`, `satellite_longitude`, `satellite_altitude` (km), `view_zenith_angle`, `view_azimuth_angle` (degrees) | (along_track, cross_track) | — | optional; if not in the files, recorded in `availability` (§6.1) |

**Dimensionality metadata** (dataset attributes): `spatial_dims =
"along_track cross_track level"`, `grid_spacing_km = [dx_along, dx_cross,
dz]` (NaN where not regular; MATLAB's `Spacing`, reordered), `grid =
"regular"` or `"native"`.

## 5. Required dataset attributes (all formats)

| Attribute | Example | Notes |
|---|---|---|
| `format`, `format_version` | `"swath"`, `1` | |
| `instrument` | `"AIRS"` | |
| `product` | `"3d"` | |
| `characterisation` | `"AIRS-3D-placeholder/0"` | id and version of the instrument's characterisation package (CP). *Agreed:* placeholders, `<INSTRUMENT>-<PRODUCT>-placeholder/0`, until the CP work sets the naming; the validator accepts them, and changing them later is compatible |
| `data` | `"raw"` / `"preprocessed"` / `"planetary-wave-filtered"` | "preprocessed" when the files were altered before loading (values interpolated, rounded ...) |
| `preprocessing` | `"daily files merged from ... granule files: temperatures rounded to whole kelvin"` | what was done to the files before loading; `""` for the provider's data as distributed |
| `heights` | `"interpolated"` / `"original"` / `"none"` | the vertical grid (`"none"`: 2-D products) |
| `time_range` | ISO times | as expanded by `time_handling` |
| `availability` | `"satellite_latitude: not_in_files; view_zenith_angle: not_in_files"` | every optional variable of the format that is not in the Dataset, and why (§6.1) |
| `source`, `atmosgw_curate_version`, `options` | | as limb |
| `files_read`, `files_not_found`, `multiple_matches` | | as limb |
| `state_counts` | `"temperature: outside_physical_range 12, missing 300"` | how many points of each variable are in each state |

## 6. Data states: absent, missing, bad, and qualified values

*Agreed:* "never existed", "missing" and "bad" are different, and Fuse and
the Parameter Engine must be able to tell them apart, along with any other
state that explains why a value is absent or what kind of value it is. They
are recorded at two levels.

### 6.1 Whole variables: `availability`

A variable defined by the format but not in the Dataset is listed in the
dataset attribute `availability`, with one of:

| State | Meaning | Example |
|---|---|---|
| `not_measured` | the instrument cannot provide it | wind from AIRS |
| `not_in_files` | it exists for this instrument, but not in the files read | AIRS satellite position (removed from the eepc-0184 copies) |
| `not_requested` | not asked for in this call | `source_file` without `file_source` |
| `not_provided` | the data supply no per-point value (the CP may describe it) | `temperature_uncertainty` |

A helper, `data.availability(ds, name)`, returns `"present"` or the state, so
users never have to parse the attribute.

### 6.2 Individual points: `<variable>_flag`

*Agreed:* each data variable has its own unsigned-integer state variable,
`<variable>_flag`, with the same dimensions, linked by the CF
`ancillary_variables` attribute and following the CF flag convention
(`flag_masks`, `flag_meanings`). Per-variable flags are needed because states
differ between variables at the same point (e.g. temperature measured,
pressure derived). Zero means a plain, good value.

The states are defined once, in a registry (`data/states.py`) shared by all
loaders: name, bit, whether the value is NaN, and a description. A new state
is a new entry; nothing else changes, and saved files stay readable because
each carries its own `flag_masks` and `flag_meanings`. uint16 gives 16
states, uint32 32 (the registry decides the type).

Initial registry (proposal), in three groups:

| Bit | State | Value | Meaning |
|---|---|---|---|
| | **No value: never there** | | |
| 1 | `no_observation` | NaN | nothing measured here: outside the swath, above or below a profile's measured values, padding between profiles, a gap between segments |
| 2 | `missing` | NaN | a measurement was expected but the file has none (fill value, missing scan) |
| 4 | `outside_requested_range` | NaN | removed because outside the height, latitude or longitude range asked for |
| | **No value: removed as bad** | | |
| 8 | `bad_geolocation` | NaN | latitude or longitude outside its range |
| 16 | `outside_physical_range` | NaN | temperature < 100 K, or > 400 K below 100 km |
| 32 | `statistical_outlier` | NaN | outside 3 x IQR of its height level (IASI screen) |
| 64 | `instrument_flag` | NaN | the instrument's own QC marks it bad |
| | **Value present, but qualified** | | |
| 128 | `derived` | value | computed, not measured or retrieved (e.g. pressure from height) |
| 256 | `filled` | value | filled from neighbours after loading (`analysis.fill.fill_missing`) |
| 512 | `extrapolated` | value | outside the source data's coverage |
| 1024 | `a_priori_dominated` | value | the retrieval has little information here (reserved: set from the CP's averaging kernels once available) |

Several bits can be set together (e.g. `derived` and `extrapolated`). With
`keep_outliers`, points in the "removed as bad" group keep their values,
and their bits are still set.

Instruments' own QC variables are passed through unchanged as their own
variables (e.g. `airs_temperature_qc`), with a `comment` naming the state
they feed.

### 6.3 Shared rules (both formats)

Applied by both loaders and by `to_profiles` (agreed 2026-11-25; checked by
`tests/test_states_conformance.py`):

- **`no_observation` vs `missing`.** A NaN inside what was measured is
  `missing`: in a profile, between its first and last value of that variable
  (a fill value, a gap); in a swath, a fill value at a scanned point. A NaN
  outside it is `no_observation`: above or below a profile's values, padding
  between profiles of different lengths, outside the swath, between segments.
- **Removed points.** A point removed as bad is removed in every data
  variable at that point, each carrying the reason (the limb loader follows
  MATLAB here; nadir derived variables take the main variable's states).
- **Qualified states** (`derived`, `extrapolated`, `filled`) are set only on
  values that exist.
- **`bad_geolocation`.** Nadir: latitude or longitude NaN or out of range.
  Limb: out of range only, as MATLAB, because some files have values at
  levels without a position (the GNSS files, pre-interpolated by C. Wright,
  lack positions at some levels; their values are kept, and the dataset says
  `data` "preprocessed").
- **Not raw data.** Files altered before loading give `data` "preprocessed"
  and a `preprocessing` text (GNSS; AIRS 3-D daily files).
- **Conversions keep provenance.** `to_profiles` carries `source_file` and
  `original_files`, and keeps the swath's `availability` reason for any
  variable it lists.

### 6.4 Where a variable comes from: `origin`

Each data variable has an `origin` attribute: `measured`, `retrieved`,
`derived` (with a `comment` saying how, e.g. "estimated from height with
atmosphere.height_to_pressure (U.S. Standard Atmosphere 1976)"), or `mixed`
when it differs between points, in which case the `derived` state marks the
points. *Agreed:* pressure in particular must say whether it is real or
created.

## 7. Uncertainty and provenance

*Agreed:* these belong in the characterisation package (CP).

- **In the CP:** static uncertainty (noise, biases, error correlations),
  spectral response, viewing geometry, provenance and independence (e.g. 2-D
  and 3-D AIRS share L1 radiances). The data link to it by the
  `characterisation` attribute.
- **In the data:** only per-point values supplied by the files, as
  `<variable>_uncertainty` with the same dimensions and units as
  `<variable>`, and a `comment` saying what kind (1-sigma random, total ...).
  If the files have none, `availability` says `not_provided`.

## 8. Conversions

- `to_profiles(swath) -> profile`: each (along_track, cross_track) point
  becomes a profile; `level` is kept; `time`, `latitude`, `longitude` are
  repeated over levels; flags, `granule`, `segment` and the swath indices are
  kept as per-profile variables (`swath_along`, `swath_cross`, 1-based). Lets
  the sampler, profile tools and Fuse treat swaths and profiles alike.
  Pressure, which the profile format requires and nadir retrievals here do
  not provide, is estimated from altitude (U.S. Standard Atmosphere 1976)
  with `origin` "derived" and the `derived` state. Only 3-D swaths convert.
  With `file_source`, `source_file` (repeated over levels), `source_row` (per
  profile) and `original_files` are kept; `source_profile` has no meaning for
  a swath and is `not_provided`; other absent variables keep the swath's
  `availability` reason.
- Others to be defined with the formats they connect.

## 9. Checking

`data.validate(ds, format)` checks dimensions, required variables, types,
units, value ranges of coordinates, required attributes, that every data
variable has a flag variable, that every NaN has a "no value" state, and that
every optional variable is either present or listed in `availability`. It
raises with a list of every problem. Every loader's tests run it on their
outputs.

## 10. Versioning

Each format's `format_version` is raised on any incompatible change (a
renamed or removed variable, changed units or dimensions). Adding optional
variables or states is compatible. Changes are listed in `CHANGELOG.md`.

## 11. Saving to disk

`data.to_netcdf(ds, path, pack=False)` writes the Dataset with zlib
compression. With `pack=True`, temperatures are stored as int16 with
`scale_factor` 0.01 K and an `add_offset` (resolution 0.01 K over a 655 K
span; NaN stored as `_FillValue`), halving the size again. That is ten times
finer than the best data (0.1 K), so it is safe for archives; float32 stays
the default because packing is lossy and gains nothing in memory. State flags
and coordinates are never packed.
