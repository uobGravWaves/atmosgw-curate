"""
Conversions between the common data formats.

``to_profiles(ds)``: a 3-D swath (format "swath") or a time-height series
(format "time_height") as profiles (format "profile"), so that tools written
for profiles (the sampler, profile operations, Fuse) can use nadir-sounder
and ground-based data too.

Inputs and outputs
------------------
``to_profiles(ds)`` -> ``xarray.Dataset`` in the profile format (its current
version; absent winds keep the input's ``availability`` reason).

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
- ``temperature``, its flag and ``temperature_uncertainty`` are kept as
  they are (copied). ``pressure``: the series' own if it has one, else
  estimated from altitude as for swaths, in the temperature's type, with
  ``origin`` "derived" and the ``derived`` state.
- ``session`` is kept per profile; ``source_profile`` and ``source_file``
  (``file_source=True``) are repeated over levels, and ``original_files`` is
  kept.
- ``availability`` and the attributes as for swaths, with ``converted_from``
  = "time_height". A Dataset without temperature (e.g. radar winds only)
  raises ``ValueError``: the profile format requires temperature.
- Not carried: ``additional_vars`` (e.g. ``T2_lidar``), which the profile
  format does not define, and winds (not converted yet: absent winds keep
  their ``availability`` reason; winds the input does hold are listed
  ``not_requested``); ``session`` is an extra per-profile variable.

Science
-------
Reshaping only, apart from the pressure estimate (``atmosphere``, U.S.
Standard Atmosphere 1976 [R1]).

Origin
------
New code (no MATLAB equivalent).

Validation
----------
``tests/test_data_convert.py``: shapes, order, values and states, file
provenance and availability, and that the result passes
``validate(ds, "profile")``.

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

"""

import numpy as np
import xarray as xr

from atmosgw.curate.atmosphere import height_to_pressure
from atmosgw.curate.data import states
from atmosgw.curate.data.formats import FORMATS, WIND_NAMES

_PER_PROFILE = ("granule", "granule_day", "segment", "source_row")


def to_profiles(ds):
    """A 3-D swath or a time_height Dataset as profiles. See the module
    docstring."""
    kind = ds.attrs.get("format")
    if kind == "swath":
        return _swath_to_profiles(ds)
    if kind == "time_height":
        return _time_height_to_profiles(ds)
    raise ValueError("to_profiles needs a Dataset in the 'swath' or "
                     "'time_height' format.")


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
    if "temperature_uncertainty" in ds:
        out["temperature_uncertainty"] = (pl, ds["temperature_uncertainty"].values
                                          .copy(), ds["temperature_uncertainty"].attrs)
    if "session" in ds:
        out["session"] = (("profile",), ds["session"].values, ds["session"].attrs)
    for name in ("source_profile", "source_file"):  # (time,) -> (profile, level)
        if name in ds:
            out[name] = (pl, over_levels(ds[name].values), ds[name].attrs)
    result = xr.Dataset(out)
    if "original_files" in ds.coords:  # (values only: not the site coordinates)
        result = result.assign_coords(original_files=(
            "file", ds["original_files"].values, ds["original_files"].attrs))

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
