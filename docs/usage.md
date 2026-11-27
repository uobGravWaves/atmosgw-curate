# Using the loaders

How to load limb-sounder and other profiles, nadir-sounder swaths,
ground-based time-height series and in-situ trajectories (radiosondes,
aircraft), and how to read what comes back. Every option is listed in the
loaders' docstrings (`help(load_limb_sounders)`, `help(load_nadir_sounders)`,
`help(load_time_height)`, `help(load_trajectories)`).

| Loader | Format, dimensions | Instruments |
|---|---|---|
| `load_limb_sounders` | profile, (profile, level) | ACE, GNSS, HIRDLS, MIPAS, MLS, SABER, SOFIE, ALOMAR; LIMS, ROMEX, Aeolus (winds); Sampler, Misc |
| `load_nadir_sounders` | swath, (along_track, cross_track, level) | AIRS (3-D, 2-D), IASI; model data in the AIRS layout |
| `load_time_height` | time_height, (time, level) | ALOMAR |
| `load_trajectories` | trajectory, (flight, sample) | DWDsondes, IAGOS, FAAM |

Install the package and set the data root first (`README.md`, "Quick
start" and "Where the data are"); the files each loader reads, how to lay
them out and where they come from are in `docs/data_supply.md`. The examples
below use:

```python
from atmosgw.curate.obs import (load_limb_sounders, load_nadir_sounders,
                                load_time_height, load_trajectories)
from atmosgw.curate.data import validate, availability, to_profiles
from atmosgw.curate.data import states
```

## Times

Every loader takes `time_range` as one time or two (start and end):
`numpy.datetime64`, ISO strings or `datetime` objects. `time_handling` says
how they are read; both ends are included.

- `"fuzzy"` (default): one time means its whole day. Two times are used as
  given, except that an end at exactly midnight includes that whole day.
- `"generous"`: whole days, from the start of the first day to the end of
  the day the end falls in.
- `"strict"`: exactly as given. One time selects only data at exactly that
  instant, which is usually nothing.

| `time_range` | `"fuzzy"` | `"generous"` | `"strict"` |
|---|---|---|---|
| `"2007-01-01"` | 1 Jan | 1 Jan | 1 Jan 00:00 only |
| `"2007-01-01T06:00"` | 1 Jan | 1 Jan | 06:00 only |
| `["2007-01-01T06:00", "2007-01-01T12:00"]` | 06:00-12:00 | 1 Jan | 06:00-12:00 |
| `["2007-01-01", "2007-01-02"]` | 1 and 2 Jan | 1 and 2 Jan | 1 Jan 00:00 to 2 Jan 00:00 |
| `["2007-01-01T06:00", "2007-01-02T18:00"]` | 1 Jan 06:00 to 2 Jan 18:00 | 1 and 2 Jan | 1 Jan 06:00 to 2 Jan 18:00 |

So `["2007-01-01", "2007-01-02"]` is two days, not one. "The end of a day"
is about 1 ms (1e-8 day) before the next midnight, as in MATLAB.

```python
ds = load_limb_sounders("2007-01-01", "SABER")                                # one day
ds = load_limb_sounders(["2007-01-01", "2007-01-03"], "SABER")                # three days
ds = load_limb_sounders(["2007-01-01T06:00", "2007-01-01T12:00"], "SABER")    # six hours
```

Profiles (limb) and granules (nadir) are selected by their mean time, so
they are never split; time-height rows are selected by their own time;
trajectories (flights) by their start time, so a flight is never split
either (a radiosonde launched at 23:45 belongs to that day).
MATLAB datenums are refused: convert them with
`atmosgw.curate.matlab_time.datenum_to_datetime`. Dates must exist:
"2022-01-34" is an error (MATLAB would carry it over to 3 February).

## Limb sounders and other profiles

Instruments (any case): `"ACE"`, `"GNSS"`, `"HIRDLS"`, `"MIPAS"`, `"MLS"`,
`"SABER"`, `"SOFIE"`, `"ALOMAR"` (lidar), `"LIMS"` (Nimbus-7, 1978-79),
`"ROMEX"` (GNSS radio occultation, September-November 2022), `"Aeolus"`
(Level 2B winds, see below), `"Sampler"` (output of the
[Wright and Hindley sampler](https://github.com/corwin365/Model-Sampling),
in MATLAB or netCDF files) and `"Misc"` (any files in the ACE MATLAB
format). LIMS, ROMEX, Aeolus and Sampler's netCDF files are new code; the
rest is ported from MATLAB `get_limbsounders`. ALOMAR here is the
MATLAB-parity path and gives a warning: for a time series, use
`load_time_height` (below).

```python
# default: interpolated to 18-60 km every 0.5 km, outliers removed
ds = load_limb_sounders("2007-01-01", "HIRDLS")
ds.temperature                     # (profile, level), K

# the file's own levels, 20-80 km only
ds = load_limb_sounders("2007-01-01", "HIRDLS", original_z=True,
                        height_range=(20, 80))

# your own height scale, a region, and where each profile came from
import numpy as np
ds = load_limb_sounders("2007-01-01", "SABER", height_scale=np.arange(20, 100.5, 1),
                        lat_range=(-30, 30), file_source=True, dx=True)

# planetary-wave-filtered data (temperature_residual, temperature_pw)
ds = load_limb_sounders("2007-01-01", "HIRDLS", hindley_pws=True)

# files of your own, named <identifier>...<yyyy>d<ddd>...mat
ds = load_limb_sounders("2007-01-01", "Misc", misc_path="/path/to/files",
                        misc_identifier="mydata")
```

Output dimensions are `profile` and `level` (and `mode` for
`temperature_pw`). `altitude` is the height scale (or the file heights with
`original_z`). `pressure` is in hPa; its `origin` attribute says whether it
was read or estimated from height. `additional_vars=["name", ...]` adds
other file variables.

**LIMS** is the one instrument accepting dates before 1979 (from
1978-10-25). It reads the previous and next day's files too (a file holds complete
orbits, so it runs into the next day; two start before their own day).
Each profile has the position of its 30 km tangent point at every level.
Geopotential height is `additional_vars=["geo_heights"]`.

**ROMEX** is a separate instrument from `"GNSS"`, which holds many of the
same occultations for those months. Profiles the files mark `is_bad` are
`instrument_flag` (kept with `instrument_qc=False`). A day takes about 20 s.

**Aeolus** gives winds, not temperature: `hlos_wind` (horizontal
line-of-sight wind, m s-1), `hlos_wind_uncertainty` and
`line_of_sight_azimuth`, product "hlos_wind" (profile format version 2 on).
Most Aeolus winds lie below 18 km, the bottom of the default height scale,
so pass your own (a warning says so):

```python
ds = load_limb_sounders("2019-01-07", "Aeolus", height_scale=np.arange(0, 25.5, 0.5))
ds = load_limb_sounders("2019-01-07", "Aeolus", aeolus_channel="mie", original_z=True)
```

`aeolus_channel` is "rayleigh" (Rayleigh-clear, default) or "mie"
(Mie-cloudy; which type is which is inferred from the files, UNVERIFIED).
Profiles are rebuilt from the files' observations, approximately (a
heuristic: `help(atmosgw.curate.obs._limb_aeolus)`); each level keeps its own time
and position. The HLOS sign is in `hlos_wind`'s `sign_convention`.

**Model data sampled as an instrument** (e.g. the SWANS DYAMOND runs):
`sample_path` points at the sampler output; `sampled_from` names the model
(else "a model (not named)"). Such output has `data` "sampled", the
attribute `sampled_from`, and `derived` values:

```python
ds = load_limb_sounders("2020-01-20", "Sampler", sampled_from="DYAMOND GEOS 3 km",
                        sample_path="/path/to/SWANS/HIRDLS/dyamond_geos3km")
```

## Nadir sounders

Instruments: `"AIRS"` (3-D retrieval; 2-D brightness temperature with
`product="2d"`) and `"IASI"` (3-D).

```python
# two hours of AIRS, regridded to a regular distance grid (135 along-track x 90
# cross-track per granule)
ds = load_nadir_sounders(["2007-01-01T00:00", "2007-01-01T02:00"], "AIRS")
ds.temperature                     # (along_track, cross_track, level), K, float32

# particular granules (day, granule number), on the instrument's own grid
ds = load_nadir_sounders(None, "AIRS", granules=[("2007-01-01", 5), ("2007-01-01", 6)],
                         original_grid=True)

# IASI, with the day/night retrieval flag
ds = load_nadir_sounders(["2010-07-01T00:00", "2010-07-01T03:00"], "IASI",
                         day_night=True)

# AIRS 2-D brightness temperature, two of the four channel groups (needs
# L. Hoffmann's 2-D granule files in the "2d" AIRS folder of datasets.toml)
ds = load_nadir_sounders(["2020-01-20T00:00", "2020-01-20T01:00"], "AIRS",
                         product="2d", channels=["4mu", "15mu_high"])
```

The swath is split into `segment`s wherever successive scans are more than
180 s apart (a missing granule, or separate times), or rows were left out
by a box. Each value's `latitude` and `longitude` are those of the native
pixels it comes from (`ds.attrs["geolocation"]`); `distance_along` and
`distance_cross` are the grid's swath coordinates, for spectral analysis. Each segment is
regridded on its own, and `distance_along` restarts at each. `granule` and
`granule_day` say where each along-track row came from. Scans with a missing
position or time (fill values in the files) form segments of their own:
with `original_grid=True` they are kept, as `bad_geolocation`; a regridded
swath leaves them out.

Model data sampled as AIRS, in the AIRS granule layout, load with
`file_path` and `sampled_from` (as for the limb loader):

```python
ds = load_nadir_sounders("2020-01-20", "AIRS", sampled_from="DYAMOND GEOS 3 km",
                         file_path="/path/to/SWANS/AIRS/dyamond_geos3km")
```

IASI has a statistical outlier screen on by default (`"rolling"`: values
more than 5 interquartile ranges beyond the quartiles of their height level
within a 10-minute window, day and night separately); AIRS has none. Choose
another with `screen=` (`"none"`, `"pooled"` as MATLAB, or a function); the
choice is recorded in the `screen` attribute.

## Ground-based time-height series

Instruments: `"ALOMAR"` (the lidar at Andoya). The data stay on their own
heights at the site: one 0.15 km grid holds every ALOMAR session, with no
interpolation.

```python
ds = load_time_height("2007-01-23", "ALOMAR")
ds.temperature                     # (time, level), K
ds.temperature_uncertainty         # as supplied by the files (kind UNVERIFIED)
ds.latitude, ds.longitude          # the site (no dimensions)
ds.session                         # which measurement session each row is from

# a few hours, 20-90 km, and where each row came from
ds = load_time_height(["2007-01-23T13:00", "2007-01-23T15:00"], "ALOMAR",
                      time_handling="strict", height_range=(20, 90),
                      file_source=True)

# the lower-resolution product and the ECMWF temperatures in the files
ds = load_time_height("2007-01-23", "ALOMAR", additional_vars=["T2_lidar", "T_ecmwf"])
```

Times can repeat: about 29% of ALOMAR rows share their time stamp with
another row, and the files do not say why. `ds.sel(time=slice(...))` works,
but aligning or merging with another Dataset on `time` does not.
`to_profiles(ds)` gives one profile per row (pressure estimated from
height). A whole archive is about 1.7 GB per variable.

## In-situ trajectories: radiosondes and aircraft

Instruments (any case): `"DWDsondes"` (the DWD high-resolution
radiosondes, 1990-2024; other sonde datasets will get names of their own),
`"IAGOS"` (in-service aircraft, IAGOS and MOZAIC, 1994-2024) and `"FAAM"`
(the FAAM research aircraft, core data). Each flight (an ascent, a flight)
is a row; its samples are in path order, each with its own position, time
and height, padded with NaN to the longest flight.

```python
ds = load_trajectories("2020-01-07", "DWDsondes")                # every station
ds = load_trajectories("2020-01-07", "dwdsondes", stations=["Greifswald", 1303])
ds.temperature                     # (flight, sample), K
ds.flight_id, ds.platform_id, ds.start_time                     # (flight)

ds = load_trajectories("2019-03-18", "IAGOS", additional_vars=["O3_P1"])
ds.temperature_source, ds.altitude_source                       # (flight): which file variables

ds = load_trajectories("2025-07-09", "FAAM")                     # 1 Hz
ds = load_trajectories("2025-07-09", "FAAM", sample_rate="full") # 32 Hz
ds.upward_wind                     # FAAM only
```

- **DWDsondes**: positions drift with the balloon; the early years' files
  give none, so the station's position is used (recorded in the
  `latitude` comment). Altitude is geopotential height. `reference_time` is
  the nominal synoptic time. `stations` takes ids or names. The files have
  no field description: meanings are inferred (`docs/data_supply.md`).
- **IAGOS**: unzipped or zipped flight files (of several copies, the
  newest version: `docs/data_supply.md`). Winds are the files'
  `zon_wind_AC`, `mer_wind_AC`. Temperature is
  the IAGOS package's where the flight has usable values (not all flagged
  erroneous or missing), else the aircraft's; altitude likewise GPS (none in
  MOZAIC-era files), else barometric. Pressure, temperature and wind flagged
  erroneous (3) or missing (7) are `instrument_flag`; "limited" and "not
  validated" are kept; position and altitude flags are not applied.
- **FAAM**: `sample_rate="1hz"` (default): each second's mean of its
  unflagged 32 Hz samples (`instrument_flag` if none); `"full"`: every
  sample, flagged ones `instrument_flag` (extras at 20 or 64 Hz need
  `"1hz"`). A flagged position has no value. The highest revision of a
  flight is read. Altitude is GPS-inertial; temperature the deiced sensor
  (`TAT_DI_R`).

Every data variable is as measured (`origin` "measured").
`help(atmosgw.curate.obs._traj_dwd_sondes)` (and `_traj_iagos`,
`_traj_faam`) say how each file field is used and what was checked.

## What comes back: formats, states and attributes

Each loader returns an `xarray.Dataset` in a named format: `"profile"` for
`load_limb_sounders`, `"swath"` for `load_nadir_sounders`, `"time_height"`
for `load_time_height`, `"trajectory"` for `load_trajectories` (defined in
`docs/design/data_formats.md`).

```python
ds.attrs["format"]                 # "profile", "swath", "time_height" or "trajectory"
validate(ds, ds.attrs["format"])   # raises FormatError listing any problem
```

**Why a value is missing.** Every data variable (`temperature`, `pressure`,
`brightness_temperature` ...) has a flag variable `<name>_flag` of the same
shape: one bit per data state, zero for a plain good value.

```python
ds = load_limb_sounders("2007-01-01", "HIRDLS", original_z=True)
flag = ds.temperature_flag.values
states.has_state(flag, "outside_physical_range")   # bool array
states.no_value_mask(flag)                         # True wherever there is no value
ds.attrs["state_counts"]                           # e.g. "temperature: no_observation 1264, ..."
```

| State | Value | Meaning |
|---|---|---|
| `no_observation` | NaN | nothing measured here: outside the swath, above or below a profile's values, padding between profiles |
| `missing` | NaN | a gap inside what was measured: a fill value at a scanned point, or between a profile's first and last value |
| `outside_requested_range` | NaN | outside the heights asked for |
| `bad_geolocation` | NaN | position out of range (nadir: also NaN, time missing, or a scan of implausible width; limb and trajectories: out of range only) |
| `outside_physical_range` | NaN | temperature < 100 K (IASI: also exactly 100 K, the retrieval's lower bound), or > 400 K (limb, time-height and trajectories: only below 100 km; nadir: at any height) |
| `statistical_outlier` | NaN | marked by the instrument's screen (IASI) |
| `derived` | value | computed, not measured or retrieved (e.g. pressure from height) |
| `extrapolated` | value | outside the data's own coverage, with `extrapolate=True` (nadir: also edge columns given a segment's edge value, Q14) |
| `instrument_flag` | NaN | failed the provider's own screening (below) |
| `a_priori_dominated` | value | mostly a priori, as the provider marks it (HIRDLS negative precision) |

Of the providers' own quality screening, these are applied as
`instrument_flag`: MLS's (its data quality documents' temperature rules),
SOFIE's endorsed heights (15-88 km), ROMEX's `is_bad`, Aeolus's `valid`
flag, IAGOS's validity flags 3 and 7, and FAAM's flags. HIRDLS values
dominated by the a priori are marked `a_priori_dominated` but kept. In
`load_limb_sounders`, `instrument_qc=False` switches all of this off;
`load_trajectories` has no such switch (IAGOS and FAAM flags are always
applied). For the other instruments, request their quality fields with
`additional_vars` and screen as the product's documentation says.

The same rules hold in every format (`docs/design/data_formats.md` section
6.3). A point removed as bad is removed in every data variable at that
point (trajectories: only in the failed sensor's variable, unless the
position is bad). The value states (`derived`, `extrapolated`) are only on
values that exist. To judge coverage, use `missing` and `no_observation`
together, not `no_observation` alone. With `keep_outliers=True` (every
loader) the values removed as bad are kept and their states still set.
The full list is `states.STATES`.

**Not raw data.** `ds.attrs["data"]` is "preprocessed" when the files were
altered before loading, and `ds.attrs["preprocessing"]` says how. This
applies to the GNSS files (pre-interpolated to a fixed height grid) and the
AIRS 3-D daily files (see "Things to know"). It is "sampled" for model data
sampled as an instrument (Sampler, or `sampled_from`).

**Variables that do not exist.** Optional variables of the format that are
not in the Dataset are listed in the `availability` attribute, with the
reason; absent temperature, pressure and winds are `not_measured` (e.g.
winds from SABER or AIRS) or `not_in_files` (each loader's docstring says
which):

```python
availability(ds, "source_file")             # "present", or e.g. "not_requested"
airs = load_nadir_sounders(["2007-01-01T00:00", "2007-01-01T01:00"], "AIRS")
availability(airs, "satellite_latitude")    # "not_in_files"
availability(airs, "eastward_wind")         # "not_measured"
```

**The record of the load**, in the attributes: `instrument`, `time_range` (as
expanded), `options` (the call's arguments, JSON), `files_read`,
`files_not_found` (warned), `files_without_data` (files that hold no
values of the main variable; warned), `multiple_matches` (each day, file,
flight or revision with more than one candidate, and the one used), `source_versions` (the data versions the files
state, with counts; with `file_source`, `original_file_version` per file),
`atmosgw_curate_version`, and per variable `source_variable` and `comment`. All are text, so they survive saving.

## Swaths and time-height series as profiles

`to_profiles` turns a 3-D swath (one profile per swath point) or a
time-height series (one profile per time) into the profile format, so
profile tools can use them too:

```python
prof = to_profiles(load_nadir_sounders(["2007-01-01T00:00", "2007-01-01T01:00"], "AIRS"))
prof.temperature                   # (profile, level)
prof.pressure                      # estimated from altitude, marked derived
prof = to_profiles(load_time_height("2007-01-23", "ALOMAR"))
```

Combining limb profiles with `to_profiles(swath)` output mixes float64 and
float32 data variables: both are valid profile data, but stacking them gives
float64.

## Saving

```python
from atmosgw.curate.data import load_netcdf, to_netcdf

to_netcdf(ds, "/path/to/output.nc")              # checked first: invalid data are refused
back = load_netcdf("/path/to/output.nc")         # the same Dataset, checked again
to_netcdf(ds, "/path/to/small.nc", pack=True)    # smaller, lossy: see below
```

`to_netcdf` compresses the file, keeps times exact (`datetime64[us]`), and
writes under a temporary name, renamed only when the file is complete, so a
failed write leaves nothing behind. An existing file is an error unless
`clobber=True`. With `pack=True`, temperatures are stored to 0.01 K,
uncertainties to 0.001 K, and computed positions and heights to about 1 cm
and 1 mm (as scaled integers); pressure, times, flags, extra variables and
positions read straight from float32 files are stored as they are. Missing
values in text extras (some GNSS `MetaData` columns) come back as empty
text. Use packing for archives, not where a result must be reproduced
bit for bit. (`ds.to_netcdf` from xarray works too, but writes invalid
Datasets as well, and times come back as `datetime64[ns]`.)

A regridded day of AIRS is 0.6 GB returned and needs about 2.6 GB while
loading (2007-01-01); IASI is similar. Save to a disk with
room (on the group's Linux machine, under `~/Data/<you>/`, not your home
folder).

## Same option, different behaviour

| Option | `load_limb_sounders` | `load_nadir_sounders` | `load_trajectories` |
|---|---|---|---|
| `lat_range`, `lon_range` | keeps profiles with a level inside the box | keeps along-track rows with any point in the box | not available |
| default heights | 18-60 km every 0.5 km (as MATLAB) | a 3 km grid from the lowest file level (as MATLAB) | none: samples at their own heights |
| `extrapolate` | vertical, beyond each profile | horizontal, in the regridding | not available |

## Things to know

- The AIRS 3-D daily files store times as float32, so times are up to 8-32 s
  off, and temperatures rounded to whole kelvin (about 0.3 K rms).
- Limb data keep the type the files give (mostly float64); nadir data are
  float32. Either is far finer than the data's precision.
- UNVERIFIED: the AIRS and IASI start dates in the loaders, and that IASI's
  retrieval switches between day and night at the same 96 degree solar zenith
  angle as AIRS.
- HIRDLS and MLS times are converted without leap seconds (a few seconds).
- MLS heights come from a standard atmosphere on MLS's pressure levels.
- GNSS: the group's merged files hold the same occultation from several
  processing centres; one copy is returned (`keep_duplicate_occultations`
  for all). JPL profiles have no per-level positions (their profile
  position is used at every level) and hour-rounded times (the occultation
  id's time is used). See Q9-Q11 in `docs/quirks.md`.
- The newer readers' data oddities and gaps (left as distributed) are in
  `docs/data_supply.md`; their open questions in their docstrings.
