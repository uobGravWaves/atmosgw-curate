"""
Conversions between the common data formats.

``to_profiles(ds)``: a 3-D swath (format "swath"), a time-height series
(format "time_height") or radiosonde ascents (format "trajectory") as
profiles (format "profile"), so that tools written for profiles (the
sampler, profile operations, Fuse, the vertical detrending filter) can use
nadir-sounder, ground-based and radiosonde data too.

Inputs and outputs
------------------
``to_profiles(ds, dz_km=0.05, height_range=None)`` -> ``xarray.Dataset`` in
the profile format (its current version; absent winds keep the input's
``availability`` reason). ``dz_km`` (float, km) and ``height_range``
((low, high) in km, or None) apply to trajectories only: given (other than
the defaults) for a swath or time_height Dataset they are a TypeError.

- Each (along_track, cross_track) point becomes a profile, in along-track
  order then cross-track (profile p = along * n_cross + cross, 0-based);
  ``level`` is kept.
- ``latitude``, ``longitude``, ``time`` are repeated over levels;
  ``altitude`` is repeated over profiles (or reshaped, if it varies by point).
- ``temperature`` and its flag (and ``temperature_uncertainty``, if
  present) are reshaped. ``pressure`` (which the profile format requires):
  the swath's own, with its flag and attributes, if it has one; otherwise
  (as for AIRS and IASI here) estimated from altitude with
  ``atmosphere.height_to_pressure`` (U.S. Standard Atmosphere 1976 [R1]):
  ``origin`` "derived", and the ``derived`` state at every point.
- Per profile: ``swath_along``, ``swath_cross`` (1-based positions in the
  swath), ``granule``, ``granule_day``, ``segment``, and ``day_retrieval`` and
  ``source_row`` if present.
- ``source_file`` (if present, from ``file_source=True``): the row's value
  repeated over cross-track positions and levels, dims (profile, level); the
  ``original_files`` coordinate is kept. ``source_profile`` has no meaning
  for a swath and is listed as ``not_provided``: ``source_row`` and
  ``swath_cross`` identify the point in its file.
- ``availability``: each tracked profile variable not in the output keeps
  the swath's reason if the swath lists it (e.g. ``not_in_files``), else
  ``not_requested``.
- Attributes are kept, with ``format``, ``format_version``,
  ``availability`` and ``state_counts`` set for the profile format, and
  ``converted_from`` = "swath".

2-D swaths have no levels and raise ``ValueError``.

From "time_height" (``load_time_height``):

- Each time (row) becomes a profile, in the same order; ``level`` is kept.
- ``latitude``, ``longitude`` (the site) and ``time`` are repeated over
  levels, ``altitude`` over profiles.
- ``temperature``, its flag and the uncertainties (``temperature_``, and
  ``pressure_uncertainty`` with the series' own pressure) are kept as they
  are (copied); references in ``ancillary_variables`` to variables not
  carried are dropped.
  ``pressure``: the series' own if it has one, else
  estimated from altitude as for swaths, in the temperature's type, with
  ``origin`` "derived" and the ``derived`` state.
- ``session`` is kept per profile (copied); ``source_profile`` and
  ``source_file`` (``file_source=True``) are repeated over levels, and
  ``original_files`` with ``original_file_version`` are kept.
- ``availability`` and the attributes as for swaths, with ``converted_from``
  = "time_height". A Dataset without temperature (e.g. radar winds only)
  raises ``ValueError``: the profile format requires temperature.
- Not carried: ``additional_vars`` (e.g. ``T2_lidar``), which the profile
  format does not define, and winds (not converted yet: absent winds keep
  their ``availability`` reason; winds the input does hold are listed
  ``not_requested``); ``session`` is an extra per-profile variable.

From "trajectory" (``load_trajectories``; radiosondes only):

- Only datasets whose ``instrument`` ends in "sondes" (any case, e.g.
  "DWDsondes"); others (IAGOS, FAAM aircraft) raise ``ValueError``. A
  Dataset without temperature raises ``ValueError``.
- Each flight becomes a profile, in the same order. Its ascent is kept: the
  samples whose altitude exceeds every earlier valid altitude (a strictly
  increasing running maximum, so it ends at the flight's highest altitude
  and drops the descent after burst and any dip below an earlier height);
  samples with NaN altitude are dropped. A flight with fewer than two such
  samples has no values.
- One regular grid for all flights: ``level`` ascending, ``altitude`` =
  k * ``dz_km`` (integer k, not repeated addition), from the lowest to the
  highest multiple of ``dz_km`` inside some flight's ascent, or the
  multiples inside ``height_range`` (its ends to 1e-9 of a step, so
  (0, 30) includes 30 km). ``altitude`` is repeated over profiles, with
  the trajectory's attributes (DWD: geopotential height).
- Interpolation, between the two ascent samples either side of each level:
  linear in height for ``temperature``, the winds, ``latitude``,
  ``longitude`` (unwrapped along the ascent where it crosses the dateline,
  wrapped back to [-180, 180] as MATLAB ``wrapTo180``) and ``time`` (offsets
  from ``start_time``, rounded to the us); ``pressure`` linearly in
  log(pressure), as the limb loader (Q6, [R3]; pressures <= 0 have no log).
  No extrapolation: outside a flight's own ascent (below its first or above
  its last sample) every value is NaN (NaT). No bridging: if either of the
  two samples of a level is NaN for a variable, that variable is NaN there.
  Each interpolated variable's ``comment`` says so. Data variables keep the
  input's type; ``latitude``, ``longitude``, ``altitude`` are float64.
- States: on values, the qualified states (and, with ``keep_outliers``,
  removal states) of the two samples; at a NaN, the removal states of the
  two samples (e.g. ``outside_physical_range``), then, as the shared rules
  (``data_formats.md`` 6.3), ``missing`` between the profile's first and
  last value of that variable and ``no_observation`` outside it (so outside
  the flight's ascent). No new qualified state is set.
- ``pressure``: the sondes' own; if absent, estimated from the grid's
  heights as for swaths (``origin`` "derived", ``derived`` state).
- Per profile: every per-flight variable (``flight_id``, ``platform_id``,
  ``platform_name``, ``platform_type``, ``start_time``, DWD
  ``reference_time``). ``source_file`` (``file_source=True``) is repeated
  over levels, and ``original_files`` with ``original_file_version`` are
  kept. ``source_profile`` has no meaning for a flight and is listed as
  ``not_provided`` (``flight_id`` and ``source_file`` identify the ascent).
- ``availability`` as for swaths (the trajectory's reasons, e.g.
  ``upward_wind`` ``not_measured``; uncertainties ``not_provided``);
  attributes kept, with ``converted_from`` = "trajectory" and
  ``conversion`` (JSON: ``dz_km``, ``height_range``, the levels and the
  rules).
- Not carried: ``additional_vars`` (e.g. ``AE_RF``), which the profile
  format does not define, and any uncertainty variable.

Science
-------
Swaths and time_height: reshaping only, apart from the pressure estimate
(``atmosphere``, U.S. Standard Atmosphere 1976 [R1]).

Trajectories: the default 50 m grid is the spacing used by Moffat-Griffin
et al. for radiosondes (the paper was not checked here:
UNVERIFIED). Pressure is interpolated in
log(pressure), exact for an isothermal layer in geopotential height (R3 of
``docs/theory/relations.md``; DWD heights are geopotential). Caveat, to
revisit: putting 1-s sonde
samples (about 5 m apart) on a 50 m grid by linear interpolation neither
averages nor filters them, so structure finer than the grid is aliased
onto it rather than removed; a layer average would be the alternative.

Origin
------
New code (no MATLAB equivalent).

Validation
----------
``tests/test_data_convert.py``: shapes, order, values and states, file
provenance and availability, and that the result passes
``validate(ds, "profile")``. ``tests/test_convert_trajectory.py``:
radiosonde ascents against known answers (linear and exponential
profiles, burst, wiggles, gaps, the dateline, the grid options).

Status
------
Proposal. AI-originated (Claude); to be reviewed by C. Wright.

Usage
-----
>>> from atmosgw.curate.data.convert import to_profiles  # doctest: +SKIP
>>> profiles = to_profiles(load_nadir_sounders(
...     ["2007-01-01T00:00", "2007-01-01T01:00"], "AIRS"))  # doctest: +SKIP

Change log
----------
2026/11/25 [Claude] initial version
2026/11/25 [Claude] keep source_file, source_row, original_files; inherit availability
2026/11/26 [Claude] time_height to profiles
2026/11/28 [Claude] review: time_height file versions, pressure uncertainty, copies
2026/11/28 [Claude] radiosonde trajectories to profiles on a regular height grid

"""

import json

import numpy as np
import xarray as xr

from atmosgw.curate.atmosphere import height_to_pressure
from atmosgw.curate.data import states
from atmosgw.curate.data.formats import FORMATS, WIND_NAMES

_PER_PROFILE = ("granule", "granule_day", "segment", "source_row")

DZ_KM = 0.05  # km: default radiosonde grid spacing (50 m; see the module docstring)
_RANGE_TOLERANCE = 1e-9  # in grid steps: a given height_range end on a multiple
_REMOVED_BITS = states.FLAG_DTYPE(sum(s.bit for s in states.STATES
                                      if s.group == "removed as bad"))
_QUALIFIED_BITS = states.FLAG_DTYPE(sum(s.bit for s in states.STATES
                                        if not s.no_value))


def to_profiles(ds, dz_km=DZ_KM, height_range=None):
    """A 3-D swath, a time_height Dataset or radiosonde trajectories as
    profiles. See the module docstring."""
    kind = ds.attrs.get("format")
    if kind == "trajectory":
        return _trajectory_to_profiles(ds, dz_km, height_range)
    if kind not in ("swath", "time_height"):
        raise ValueError("to_profiles needs a Dataset in the 'swath', "
                         "'time_height' or 'trajectory' format.")
    # (options of trajectories only: a value other than the default, as the
    # loaders treat misplaced options)
    for option, given in (("dz_km", dz_km != DZ_KM),
                          ("height_range", height_range is not None)):
        if given:
            raise TypeError(f"{option} applies to trajectories (radiosondes) only, "
                            f"not to a {kind} Dataset.")
    if kind == "swath":
        return _swath_to_profiles(ds)
    return _time_height_to_profiles(ds)


def _swath_to_profiles(ds):
    """A 3-D swath as profiles."""
    if "level" not in ds.dims or "temperature" not in ds:
        raise ValueError("to_profiles needs a 3-D swath (with levels).")
    n_along, n_cross = ds.sizes["along_track"], ds.sizes["cross_track"]
    n_lev = ds.sizes["level"]
    n_prof = n_along * n_cross
    pl = ("profile", "level")

    def per_point(name):  # (along, cross) -> (profile,)
        return ds[name].values.reshape(n_prof)

    def over_levels(values):  # (profile,) -> (profile, level)
        return np.repeat(values[:, None], n_lev, axis=1)

    altitude = ds["altitude"].values
    if altitude.ndim == 1:
        altitude = np.broadcast_to(altitude, (n_prof, n_lev)).copy()
    else:
        altitude = altitude.reshape(n_prof, n_lev)
    temp = ds["temperature"].values.reshape(n_prof, n_lev).copy()
    temp_flag = ds["temperature_flag"].values.reshape(n_prof, n_lev).copy()
    pressure, pressure_flag, pressure_attrs = _pressure(
        ds, altitude, np.float32, lambda name: ds[name].values.reshape(n_prof, n_lev))

    along, cross = np.meshgrid(np.arange(1, n_along + 1), np.arange(1, n_cross + 1),
                               indexing="ij")
    out = {
        "latitude": (pl, over_levels(per_point("latitude")), ds["latitude"].attrs),
        "longitude": (pl, over_levels(per_point("longitude")), ds["longitude"].attrs),
        "time": (pl, over_levels(per_point("time")), ds["time"].attrs),
        "altitude": (pl, altitude.astype(np.float64), ds["altitude"].attrs),
        "temperature": (pl, temp, ds["temperature"].attrs),
        "temperature_flag": (pl, temp_flag, ds["temperature_flag"].attrs),
        "pressure": (pl, pressure, pressure_attrs),
        "pressure_flag": (pl, pressure_flag, states.flag_attrs()),
        "swath_along": (("profile",), along.reshape(n_prof).astype(np.int32),
                        {"long_name": "along-track position in the swath (1-based)"}),
        "swath_cross": (("profile",), cross.reshape(n_prof).astype(np.int32),
                        {"long_name": "cross-track position in the swath (1-based)"}),
    }
    for name in _PER_PROFILE:
        if name in ds:
            out[name] = (("profile",), np.repeat(ds[name].values, n_cross),
                         ds[name].attrs)
    if "day_retrieval" in ds:
        out["day_retrieval"] = (("profile",), per_point("day_retrieval"),
                                ds["day_retrieval"].attrs)
    if "temperature_uncertainty" in ds:  # (supplied: kept)
        out["temperature_uncertainty"] = (pl, ds["temperature_uncertainty"].values
                                          .reshape(n_prof, n_lev).copy(),
                                          ds["temperature_uncertainty"].attrs)
    if "source_file" in ds:  # row value over cross-track positions and levels
        out["source_file"] = (pl, over_levels(np.repeat(ds["source_file"].values,
                                                        n_cross)),
                              ds["source_file"].attrs)
    result = xr.Dataset(out)
    if "original_files" in ds.coords:
        result = result.assign_coords(original_files=ds["original_files"])

    # absent tracked variables keep the swath's reason; source_profile has no
    # meaning for a swath (source_row and swath_cross locate the point)
    swath_availability = states.parse_availability(ds.attrs.get("availability", ""))
    tracked = [n for n, s in FORMATS["profile"]["variables"].items()
               if s.get("tracked") and n not in result]
    availability = {n: swath_availability.get(n, "not_requested") for n in tracked}
    for name in WIND_NAMES:  # (a swath saved before format version 2 lists none)
        if name in availability and name not in swath_availability:
            availability[name] = "not_measured"
    if "source_profile" in availability:
        availability["source_profile"] = "not_provided"
    for name in ("temperature_uncertainty", "pressure_uncertainty"):
        if name not in result:
            availability[name] = swath_availability.get(name, "not_provided")
    result.attrs = dict(ds.attrs)
    for name in ("spatial_dims", "grid_spacing_km", "grid"):
        result.attrs.pop(name, None)
    result.attrs.update({
        "format": "profile", "format_version": FORMATS["profile"]["version"],
        "converted_from": "swath",
        "availability": states.format_availability(availability),
        "state_counts": states.counts_text(
            {n: result[f"{n}_flag"].values for n in ("temperature", "pressure")}),
    })
    return result


def _time_height_to_profiles(ds):
    """A time_height Dataset (with temperature) as profiles: one per time."""
    if "temperature" not in ds:
        raise ValueError("to_profiles needs temperature (the profile format "
                         "requires it).")
    n_prof, n_lev = ds.sizes["time"], ds.sizes["level"]
    pl = ("profile", "level")

    def over_levels(values):  # (profile,) -> (profile, level)
        return np.repeat(np.asarray(values)[:, None], n_lev, axis=1)

    altitude = np.broadcast_to(ds["altitude"].values, (n_prof, n_lev)).astype(
        np.float64)
    temp = ds["temperature"].values.copy()  # (a copy: the input is not shared)
    pressure, pressure_flag, pressure_attrs = _pressure(
        ds, altitude, temp.dtype, lambda name: ds[name].values.copy())

    def site(name):  # the site position at every point
        return np.full((n_prof, n_lev), float(ds[name].values))

    out = {
        "latitude": (pl, site("latitude"), ds["latitude"].attrs),
        "longitude": (pl, site("longitude"), ds["longitude"].attrs),
        "time": (pl, over_levels(ds["time"].values), ds["time"].attrs),
        "altitude": (pl, altitude, ds["altitude"].attrs),
        "temperature": (pl, temp, ds["temperature"].attrs),
        "temperature_flag": (pl, ds["temperature_flag"].values.copy(),
                             ds["temperature_flag"].attrs),
        "pressure": (pl, pressure, pressure_attrs),
        "pressure_flag": (pl, pressure_flag, states.flag_attrs()),
    }
    own_pressure = "pressure" in ds and "pressure_flag" in ds  # (as _pressure)
    carried = ("temperature_uncertainty",) + (("pressure_uncertainty",)
                                              if own_pressure else ())
    for name in carried:
        if name in ds:
            out[name] = (pl, ds[name].values.copy(), ds[name].attrs)
    if "session" in ds:  # (a copy: the input is not shared)
        out["session"] = (("profile",), ds["session"].values.copy(),
                          ds["session"].attrs)
    for name in ("source_profile", "source_file"):  # (time,) -> (profile, level)
        if name in ds:
            out[name] = (pl, over_levels(ds[name].values), ds[name].attrs)
    result = xr.Dataset(out)
    # the files and the version each states (values only: not the site
    # coordinates)
    files = {name: ("file", ds[name].values.copy(), ds[name].attrs)
             for name in ("original_files", "original_file_version") if name in ds}
    if files:
        result = result.assign_coords(files)
    # ancillary references only to variables carried over
    for name in ("temperature", "pressure"):
        refs = str(result[name].attrs.get("ancillary_variables", "")).split()
        if refs:
            kept = [r for r in refs if r in result]
            attrs = {k: v for k, v in result[name].attrs.items()
                     if k != "ancillary_variables"}
            if kept:
                attrs["ancillary_variables"] = " ".join(kept)
            result[name].attrs = attrs

    # absent tracked variables keep the time_height reason where it has one
    th_availability = states.parse_availability(ds.attrs.get("availability", ""))
    tracked = [n for n, s in FORMATS["profile"]["variables"].items()
               if s.get("tracked") and n not in result]
    availability = {n: th_availability.get(n, "not_requested") for n in tracked}
    if "pressure_uncertainty" not in result:
        availability["pressure_uncertainty"] = "not_provided"
    if "temperature_uncertainty" not in result:
        availability["temperature_uncertainty"] = th_availability.get(
            "temperature_uncertainty", "not_provided")
    result.attrs = dict(ds.attrs)
    result.attrs.update({
        "format": "profile", "format_version": FORMATS["profile"]["version"],
        "converted_from": "time_height",
        "availability": states.format_availability(availability),
        "state_counts": states.counts_text(
            {n: result[f"{n}_flag"].values for n in ("temperature", "pressure")}),
    })
    return result


def _trajectory_to_profiles(ds, dz_km, height_range):
    """Radiosonde ascents (trajectory format) as profiles on one regular
    height grid: one profile per flight."""
    instrument = str(ds.attrs.get("instrument", ""))
    if not instrument.lower().endswith("sondes"):
        raise ValueError(f"to_profiles can put only radiosonde ascents on a height "
                         f"grid; instrument {instrument!r} is not a radiosonde "
                         f"dataset (aircraft flights are not converted).")
    if "temperature" not in ds:
        raise ValueError("to_profiles needs temperature (the profile format "
                         "requires it).")
    dz_km, height_range = _check_grid_options(dz_km, height_range)
    alt = np.asarray(ds["altitude"].values, dtype=np.float64)
    keep = _ascent(alt)
    grid = _height_grid(alt, keep, dz_km, height_range)
    lower, upper, weight, inside = _interpolation_plan(alt, keep, grid)
    n_prof, n_lev = alt.shape[0], grid.size
    pl = ("profile", "level")

    def take(values, index):  # (flight, sample) -> (profile, level)
        values = np.asarray(values)
        if values.shape[1] == 0:
            return np.zeros(index.shape, dtype=values.dtype)
        return np.take_along_axis(values, index, axis=1)

    def interpolate(values):  # linear in height; NaN outside the ascent's range
        values = np.asarray(values, dtype=np.float64)
        a, b = take(values, lower), take(values, upper)
        out = a + weight * (b - a)  # (NaN if either neighbour is NaN)
        out[~inside] = np.nan
        return out

    def grid_flag(name, value):
        source = ds[f"{name}_flag"].values
        carried = np.where(inside, take(source, lower) | take(source, upper), 0)
        return _grid_flags(value, carried.astype(states.FLAG_DTYPE))

    how = (f"linearly interpolated in height{{}} from the ascent's samples to a "
           f"regular {dz_km:g} km grid by data.convert.to_profiles (no averaging)")

    def noted(attrs, extra=""):  # the variable's attributes, with the note added
        attrs = dict(attrs)
        prior = attrs.get("comment")
        text = how.format(extra)
        attrs["comment"] = f"{prior}; {text}" if prior else text
        return attrs

    altitude = np.broadcast_to(grid, (n_prof, n_lev)).astype(np.float64)
    # longitude: unwrapped along each ascent where it crosses the dateline,
    # interpolated, and wrapped back to [-180, 180]
    lon = np.array(ds["longitude"].values, dtype=np.float64)
    for k in range(n_prof):
        valid = keep[k] & np.isfinite(lon[k])
        if valid.sum() > 1 and np.any(np.abs(np.diff(lon[k, valid])) > 180):
            lon[k, valid] = np.unwrap(lon[k, valid], period=360.0)
    out = {
        "latitude": (pl, interpolate(ds["latitude"].values),
                     noted(ds["latitude"].attrs)),
        "longitude": (pl, _wrap_to_180(interpolate(lon)),
                      noted(ds["longitude"].attrs)),
        "time": (pl, _interpolate_time(ds, interpolate, (n_prof, n_lev)),
                 noted(ds["time"].attrs)),
        "altitude": (pl, altitude, ds["altitude"].attrs),
    }
    data_names = [n for n in ("temperature", "pressure") + WIND_NAMES
                  if n in ds and f"{n}_flag" in ds]
    gridded = {}
    for name in data_names:
        dtype = ds[name].dtype
        if name == "pressure":  # in log(pressure), as the limb loader (Q6, [R3])
            with np.errstate(invalid="ignore", divide="ignore"):
                log_p = np.log(np.where(ds[name].values > 0, ds[name].values, np.nan))
            value = np.exp(interpolate(log_p))
            attrs = noted(ds[name].attrs, " (pressure: in log(pressure))")
        else:
            value = interpolate(ds[name].values)
            attrs = noted(ds[name].attrs)
        gridded[name] = value.astype(dtype)
        gridded[f"{name}_flag"] = grid_flag(name, value)
        out[name] = (pl, gridded[name], attrs)
        out[f"{name}_flag"] = (pl, gridded[f"{name}_flag"], states.flag_attrs())
    if "pressure" not in gridded:  # (estimated from the grid's heights)
        pressure, pressure_flag, pressure_attrs = _pressure(
            ds, altitude, ds["temperature"].dtype, lambda name: gridded[name])
        out["pressure"] = (pl, pressure, pressure_attrs)
        out["pressure_flag"] = (pl, pressure_flag, states.flag_attrs())
        data_names.insert(1, "pressure")
    for name in ds.variables:  # per flight -> per profile (not source_file)
        if ds[name].dims == ("flight",) and name != "source_file":
            out[name] = (("profile",), ds[name].values.copy(), ds[name].attrs)
    if "source_file" in ds:  # (flight,) -> (profile, level)
        out["source_file"] = (pl, np.repeat(ds["source_file"].values[:, None], n_lev,
                                            axis=1), ds["source_file"].attrs)
    result = xr.Dataset(out)
    files = {name: ("file", ds[name].values.copy(), ds[name].attrs)
             for name in ("original_files", "original_file_version") if name in ds}
    if files:
        result = result.assign_coords(files)

    # absent tracked variables keep the trajectory's reason where it has one;
    # source_profile has no meaning for a flight (flight_id and source_file
    # identify it); uncertainties as the trajectory lists them
    traj_availability = states.parse_availability(ds.attrs.get("availability", ""))
    tracked = [n for n, s in FORMATS["profile"]["variables"].items()
               if s.get("tracked") and n not in result]
    availability = {n: traj_availability.get(n, "not_requested") for n in tracked}
    if "source_profile" in availability:
        availability["source_profile"] = "not_provided"
    for name in data_names:
        uncertainty = f"{name}_uncertainty"
        if uncertainty not in result:
            availability[uncertainty] = traj_availability.get(uncertainty,
                                                              "not_provided")
    result.attrs = dict(ds.attrs)
    result.attrs.update({
        "format": "profile", "format_version": FORMATS["profile"]["version"],
        "converted_from": "trajectory",
        "conversion": json.dumps({
            "dz_km": dz_km,
            "height_range": None if height_range is None else list(height_range),
            "levels": int(n_lev),
            "lowest_km": float(grid[0]) if n_lev else None,
            "highest_km": float(grid[-1]) if n_lev else None,
            "ascent": "samples above every earlier altitude (to the highest)",
            "interpolation": "linear in height; pressure linear in log(pressure); "
                             "no extrapolation; not across a missing sample"}),
        "availability": states.format_availability(availability),
        "state_counts": states.counts_text(
            {n: result[f"{n}_flag"].values for n in data_names}),
    })
    return result


def _check_grid_options(dz_km, height_range):
    """(dz_km, height_range) as floats, checked."""
    try:
        dz_km = float(dz_km)
    except (TypeError, ValueError):
        raise ValueError(f"dz_km must be a number (km), not {dz_km!r}.") from None
    if not np.isfinite(dz_km) or dz_km <= 0:
        raise ValueError(f"dz_km must be a positive number of km, not {dz_km!r}.")
    if height_range is None:
        return dz_km, None
    try:
        low, high = (float(h) for h in height_range)
    except (TypeError, ValueError):
        raise ValueError("height_range must be (low, high) in km.") from None
    if not (np.isfinite(low) and np.isfinite(high) and low < high):
        raise ValueError(f"height_range must be (low, high) in km with low < high, "
                         f"not {height_range!r}.")
    return dz_km, (low, high)


def _ascent(alt):
    """[flight x sample] bool: the ascent of each flight, the samples whose
    altitude exceeds every earlier valid altitude (a strictly increasing
    running maximum, so ending at the highest); NaN altitudes are dropped."""
    valid = np.isfinite(alt)
    running = np.maximum.accumulate(np.where(valid, alt, -np.inf), axis=1)
    earlier = np.concatenate([np.full((alt.shape[0], 1), -np.inf),
                              running[:, :-1]], axis=1)[:, :alt.shape[1]]
    return valid & (alt > earlier)


def _height_grid(alt, keep, dz_km, height_range):
    """The levels k * dz_km (km, ascending, integer k): those inside
    height_range (to _RANGE_TOLERANCE of a step at its ends), else from the
    lowest to the highest multiple of dz_km inside some flight's ascent."""
    if height_range is not None:
        k_low = int(np.ceil(height_range[0] / dz_km - _RANGE_TOLERANCE))
        k_high = int(np.floor(height_range[1] / dz_km + _RANGE_TOLERANCE))
        return np.arange(k_low, k_high + 1) * dz_km
    usable = keep.sum(axis=1) > 1  # (two samples are needed to interpolate)
    if not usable.any():
        return np.zeros(0)
    ascents = np.where(keep, alt, np.nan)[usable]
    low, high = np.nanmin(ascents), np.nanmax(ascents)
    k_low, k_high = int(np.ceil(low / dz_km)), int(np.floor(high / dz_km))
    # (the product k * dz_km, not the quotient, must lie inside the ascents)
    while k_low * dz_km < low:
        k_low += 1
    while (k_low - 1) * dz_km >= low:
        k_low -= 1
    while k_high * dz_km > high:
        k_high -= 1
    while (k_high + 1) * dz_km <= high:
        k_high += 1
    return np.arange(k_low, k_high + 1) * dz_km


def _interpolation_plan(alt, keep, grid):
    """(lower, upper, weight, inside), each [flight x level]: the samples
    either side of each level in the flight's ascent (indices along sample),
    the weight of the upper one, and whether the level is inside the ascent's
    height range (first to last sample, ends included)."""
    shape = (alt.shape[0], grid.size)
    lower, upper = np.zeros(shape, dtype=np.intp), np.zeros(shape, dtype=np.intp)
    weight, inside = np.zeros(shape), np.zeros(shape, dtype=bool)
    for k in range(alt.shape[0]):
        index = np.flatnonzero(keep[k])
        if index.size < 2:
            continue
        z = alt[k, index]  # (strictly increasing)
        within = (grid >= z[0]) & (grid <= z[-1])
        # segment j holds z[j] <= level < z[j+1]; the top sample goes in the last
        j = np.clip(np.searchsorted(z, grid, side="right") - 1, 0, z.size - 2)
        lower[k], upper[k] = index[j], index[j + 1]
        weight[k] = np.where(within, (grid - z[j]) / (z[j + 1] - z[j]), 0.0)
        inside[k] = within
    return lower, upper, weight, inside


def _interpolate_time(ds, interpolate, shape):
    """time on the grid (datetime64[us], NaT where there is none): offsets
    from each flight's start_time (us, float64: exact for flights under
    285 years) interpolated and rounded to the us."""
    time = ds["time"].values.astype("datetime64[us]")
    ref = ds["start_time"].values.astype("datetime64[us]")
    ref = np.where(np.isnat(ref), np.datetime64(0, "us"), ref)
    offset = interpolate((time - ref[:, None]) / np.timedelta64(1, "us"))
    out = np.full(shape, np.datetime64("NaT", "us"))
    finite = np.isfinite(offset)
    ref_grid = np.broadcast_to(ref[:, None], shape)
    out[finite] = ref_grid[finite] + np.round(offset[finite]).astype(np.int64).astype(
        "timedelta64[us]")
    return out


def _wrap_to_180(lon):
    """Longitudes (degrees) wrapped to [-180, 180] as MATLAB wrapTo180:
    values inside unchanged, positive odd multiples of 180 to 180."""
    lon = np.array(lon, dtype=np.float64)
    with np.errstate(invalid="ignore"):
        outside = (lon < -180) | (lon > 180)
    shifted = lon[outside] + 180
    wrapped = np.mod(shifted, 360)
    wrapped[(wrapped == 0) & (shifted > 0)] = 360
    lon[outside] = wrapped - 180
    return lon


def _grid_flags(value, carried):
    """The flag of a gridded variable: on values, the qualified and removal
    states of the two samples either side (carried); at a NaN, their removal
    states, then missing between the profile's first and last value along
    level (removed points count as measured) and no_observation outside."""
    nan = np.isnan(value)
    flag = np.where(nan, carried & _REMOVED_BITS,
                    carried & (_QUALIFIED_BITS | _REMOVED_BITS)).astype(
                        states.FLAG_DTYPE)
    measured = ~nan | ((flag & _REMOVED_BITS) != 0)
    cover = np.zeros(measured.shape, dtype=bool)
    if measured.shape[1]:
        levels = np.arange(measured.shape[1])
        first = np.argmax(measured, axis=1)
        last = measured.shape[1] - 1 - np.argmax(measured[:, ::-1], axis=1)
        cover = (measured.any(axis=1)[:, None] & (levels >= first[:, None])
                 & (levels <= last[:, None]))
    unexplained = nan & ~states.no_value_mask(flag)
    states.set_state(flag, unexplained & cover, "missing")
    states.set_state(flag, unexplained & ~cover, "no_observation")
    return flag


def _pressure(ds, altitude, dtype, take):
    """(pressure, its flag, its attributes) as profiles: the Dataset's own if
    it has one (take(name) reshapes a variable), else estimated from altitude
    with the U.S. Standard Atmosphere 1976, marked derived."""
    if "pressure" in ds and "pressure_flag" in ds:
        return take("pressure"), take("pressure_flag"), dict(ds["pressure"].attrs)
    pressure = np.asarray(height_to_pressure(altitude), dtype=dtype)
    flag = states.new_flags(pressure.shape)
    states.set_state(flag, np.isfinite(pressure), "derived")
    states.set_state(flag, ~np.isfinite(pressure), "no_observation")
    return pressure, flag, {
        "units": "hPa", "standard_name": "air_pressure", "origin": "derived",
        "ancillary_variables": "pressure_flag",
        "comment": "estimated from altitude with atmosphere.height_to_pressure "
                   "(U.S. Standard Atmosphere 1976) by data.convert.to_profiles"}
