# Design: common data formats

**Status:** draft for agreement, 2026-11-25. Nothing here is implemented yet.

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
as data), `flags.py` (the quality-flag registry) and `validate(ds, format)`.

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

- Times: `datetime64[us]`, UTC, named `time`. Missing: NaT.
- Latitude: degrees north, [-90, 90]. Longitude: degrees east, [-180, 180].
- Altitude: km (geometric, as given by each instrument's data; a variable
  `comment` says if it is geopotential). Pressure: hPa. Temperature: K.
- Missing values: NaN (NaT for times). No fill values in outputs.
- Source indices (rows, profiles, files): 1-based, as stored in the data.
- Signs and wave conventions: as `CONTRIBUTING.md` (upward propagation
  m <= 0; wavenumbers as 1/wavelength in km^-1; bearings clockwise from
  north).
- Names: lower case, words separated by `_`, spelled out (`temperature`, not
  `T`).
- Attributes are text, numbers or lists of them, so a Dataset survives
  `to_netcdf` unchanged.

## 3. Format "profile", version 1

**Dimensions:** `profile`, `level`; optional `mode` (planetary-wave modes),
`file`. `profile` may have length 0.

| Variable | Dimensions | Units | Required |
|---|---|---|---|
| `latitude`, `longitude` | (profile, level) | degrees | yes |
| `time` | (profile, level) | — | yes |
| `altitude` | (profile, level) | km | yes |
| `temperature` | (profile, level) | K | yes, unless PW data |
| `pressure` | (profile, level) | hPa | yes (estimated from height if absent; `comment` says so) |
| `temperature_residual`, `temperature_pw` | (profile, level[, mode]) | K | PW data |
| `quality_flag` | (profile, level) | — | yes (from version 1 of this spec; see §6) |
| `source_profile`, `source_file` | (profile, level) | — | optional |
| `original_files` | (file) | — | optional |
| `dx` | (profile, level) | km | optional |

This is the current `load_limb_sounders` output, plus `quality_flag` and the
attributes of §5. Adding these to the limb loader is a separate change; its
values do not change.

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
| `quality_flag` | as the data variable | — | yes |
| `granule`, `granule_day` | (along_track) | — | yes: the granule each row came from, and its day |
| `segment` | (along_track) | — | yes: 1-based number of the unbroken stretch the row is in |
| `distance_along`, `distance_cross` | (along_track), (cross_track) | km | yes on a regular grid |
| `source_file`, `source_row` | (along_track) | — | optional |
| `original_files` | (file) | — | optional |
| `day_retrieval` | (along_track, cross_track) | — | optional |
| native geolocation: `native_latitude`, `native_longitude`, `native_time`, `scan_position` | (native_along, native_cross) | — | optional |
| viewing geometry: `satellite_latitude`, `satellite_longitude`, `satellite_altitude` (km), `view_zenith_angle`, `view_azimuth_angle` (degrees) | (along_track, cross_track) | — | optional, when the files hold them |

Data variables are float32; coordinates float64.

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
| `characterisation` | `"AIRS-3D-Hoffmann/1"` | id and version of the instrument's characterisation package (CP) |
| `data` | `"raw"` / `"planetary-wave-filtered"` | as limb |
| `time_range` | ISO times | as expanded by `time_handling` |
| `source`, `atmosgw_curate_version`, `options` | | as limb |
| `files_read`, `files_not_found`, `multiple_matches` | | as limb |
| `flag_counts` | `"outside_physical_range: 12; filled: 12"` | how many points carry each flag |

## 6. Quality flags

One integer variable, `quality_flag`, with the same dimensions as the data,
following the CF flag convention: attributes `flag_masks` (the bit values)
and `flag_meanings` (their names, space-separated). Zero means no flag set.

The bits are defined once, in a registry (`data/flags.py`), shared by all
loaders: name, bit, and description. A new flag is a new entry; nothing else
changes, and old files remain readable because each file carries its own
`flag_masks` and `flag_meanings`. A 32-bit integer gives 32 flags; 64-bit if
ever needed.

Initial registry (proposal):

| Bit | Name | Meaning | Set by |
|---|---|---|---|
| 1 | `bad_geolocation` | latitude or longitude outside its range | both loaders (outliers) |
| 2 | `outside_physical_range` | temperature < 100 K, or > 400 K below 100 km | both loaders |
| 4 | `outside_height_scale` | height outside the requested height scale | limb |
| 8 | `statistical_outlier` | outside 3 x IQR of its height level | IASI screen |
| 16 | `filled` | value replaced by a neighbour or interpolation | IASI screen |
| 32 | `extrapolated` | outside the source data's coverage | regridding, vertical interpolation |
| 64 | `instrument_flag` | the instrument's own QC marks it bad | readers, where files have QC |

Instruments' own QC variables are passed through unchanged as their own
variables (e.g. `airs_temperature_qc`), with a `comment` naming the bit they
feed.

Unless `keep_outliers` is set, points with any of bits 1, 2, 4 or 8 are NaN,
as now; the flag says why. `filled` points keep their (filled) value, as
MATLAB.

## 7. Uncertainty and provenance

- **In the characterisation package (CP):** static uncertainty (noise,
  biases, error correlations), spectral response, viewing geometry,
  provenance and independence (e.g. 2-D and 3-D AIRS share L1 radiances).
  The data link to it by the `characterisation` attribute.
- **In the data:** only per-point values supplied by the files, as
  `<variable>_uncertainty` with the same dimensions and units as
  `<variable>`, and a `comment` saying what kind (1-sigma random, total ...).
  Absent when the files have none.

## 8. Conversions

- `to_profiles(swath) -> profile`: each (along_track, cross_track) point
  becomes a profile; `level` is kept; `time`, `latitude`, `longitude` are
  repeated over levels; `granule`, `segment` and the swath indices are kept
  as per-profile variables (`swath_along`, `swath_cross`, 1-based). Lets the
  sampler, profile tools and Fuse treat swaths and profiles alike.
- Others to be defined with the formats they connect.

## 9. Checking

`data.validate(ds, format)` checks dimensions, required variables, types,
units, value ranges of coordinates, required attributes and the flag
attributes, and raises with a list of every problem. Every loader's tests run
it on their outputs.

## 10. Versioning

Each format's `format_version` is raised on any incompatible change (a
renamed or removed variable, changed units or dimensions). Adding optional
variables or flag bits is compatible. Changes are listed in `CHANGELOG.md`.
