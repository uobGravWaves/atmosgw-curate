"""
Conversions between the common data formats.

Purpose
-------
``to_profiles(ds)``: a 3-D swath (format "swath") as profiles (format
"profile"), so that tools written for profiles (the sampler, profile
operations, Fuse) can use nadir-sounder data too.

Inputs and outputs
------------------
``to_profiles(ds)`` -> ``xarray.Dataset`` in the profile format, version 1.

- Each (along_track, cross_track) point becomes a profile, in along-track
  order then cross-track (profile p = along * n_cross + cross, 0-based);
  ``level`` is kept.
- ``latitude``, ``longitude``, ``time`` are repeated over levels;
  ``altitude`` is repeated over profiles (or reshaped, if it varies by point).
- ``temperature`` and its flag are reshaped. ``pressure`` (which the profile
  format requires and nadir retrievals here do not provide) is estimated from
  altitude with ``atmosphere.height_to_pressure`` (U.S. Standard Atmosphere
  1976 [R1]): ``origin`` "derived", and the ``derived`` state at every point.
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

Science
-------
No scientific processing beyond reshaping, except the pressure estimate
(see ``atmosphere``).

Origin
------
New code (no MATLAB equivalent).

Validation
----------
``tests/test_data_convert.py``: shapes, order, values and states, file
provenance and availability, and that the result passes
``validate(ds, "profile")``.

AI involvement
--------------
AI-originated (Claude); to be reviewed by C. Wright.

Usage
-----
>>> from atmosgw.curate.data.convert import to_profiles  # doctest: +SKIP
>>> profiles = to_profiles(load_nadir_sounders("2007-01-01", "AIRS"))  # doctest: +SKIP

Status
------
Proposal.

Change log
----------
2026/11/25 [Claude] initial version
2026/11/25 [Claude] keep source_file, source_row, original_files; inherit availability

How to cite
-----------
Cite the package (see ``CITATION.cff``).
"""

import numpy as np
import xarray as xr

from atmosgw.curate.atmosphere import height_to_pressure
from atmosgw.curate.data import states
from atmosgw.curate.data.formats import FORMATS

_PER_PROFILE = ("granule", "granule_day", "segment", "source_row")


def to_profiles(ds):
    """A 3-D swath as profiles. See the module docstring."""
    if ds.attrs.get("format") != "swath":
        raise ValueError("to_profiles needs a Dataset in the 'swath' format.")
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
    temp = ds["temperature"].values.reshape(n_prof, n_lev)
    temp_flag = ds["temperature_flag"].values.reshape(n_prof, n_lev)
    pressure = np.asarray(height_to_pressure(altitude), dtype=np.float32)
    pressure_flag = states.new_flags(pressure.shape)
    states.set_state(pressure_flag, np.isfinite(pressure), "derived")
    states.set_state(pressure_flag, ~np.isfinite(pressure), "no_observation")

    along, cross = np.meshgrid(np.arange(1, n_along + 1), np.arange(1, n_cross + 1),
                               indexing="ij")
    out = {
        "latitude": (pl, over_levels(per_point("latitude")), ds["latitude"].attrs),
        "longitude": (pl, over_levels(per_point("longitude")), ds["longitude"].attrs),
        "time": (pl, over_levels(per_point("time")), ds["time"].attrs),
        "altitude": (pl, altitude.astype(np.float64), ds["altitude"].attrs),
        "temperature": (pl, temp, ds["temperature"].attrs),
        "temperature_flag": (pl, temp_flag, ds["temperature_flag"].attrs),
        "pressure": (pl, pressure, {
            "units": "hPa", "standard_name": "air_pressure", "origin": "derived",
            "ancillary_variables": "pressure_flag",
            "comment": "estimated from altitude with atmosphere.height_to_pressure "
                       "(U.S. Standard Atmosphere 1976) by data.convert.to_profiles"}),
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
    if "source_profile" in availability:
        availability["source_profile"] = "not_provided"
    availability.update({"temperature_uncertainty": "not_provided",
                         "pressure_uncertainty": "not_provided"})
    counts = []
    for name in ("temperature", "pressure"):
        c = states.state_counts(result[f"{name}_flag"].values)
        counts.append(f"{name}: " + (", ".join(f"{k} {n}" for k, n in c.items())
                                     or "none"))
    result.attrs = dict(ds.attrs)
    for name in ("spatial_dims", "grid_spacing_km", "grid"):
        result.attrs.pop(name, None)
    result.attrs.update({
        "format": "profile", "format_version": FORMATS["profile"]["version"],
        "converted_from": "swath",
        "availability": states.format_availability(availability),
        "state_counts": "; ".join(counts),
    })
    return result
