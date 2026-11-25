"""
The common data formats, and a checker for them.

Purpose
-------
Every loader returns an ``xarray.Dataset`` in a named, versioned format, so
that the Detrender, Stockwell wrappers, Parameter Engine, sampler and Fuse can
use any instrument's data the same way. This module holds the format
definitions as data, a checker (``validate``) and the ``availability`` helper.
The formats are described in ``docs/design/data_formats.md``.

Inputs and outputs
------------------
============================================  ==================================
Function                                      Returns
============================================  ==================================
``validate(ds, format_name)``                 None; raises ``FormatError``
                                              listing every problem
``problems(ds, format_name)``                 list of str (empty if valid)
``availability(ds, name)``                    ``"present"``, or the variable's
                                              availability state
``data_variables(format_name)``               names of the format's data
                                              variables
============================================  ==================================

Formats: ``"swath"`` and ``"profile"``, version 1 (``FORMATS``).

What ``validate`` checks:

- ``format`` and ``format_version`` attributes;
- dimensions, and the order of each variable's dimensions;
- required variables (some depend on ``product``, ``grid`` or the data:
  e.g. ``temperature`` for 3-D swaths, ``brightness_temperature`` for 2-D);
- types: data variables float32; coordinates float64; ``time``
  ``datetime64[us]``; flags ``states.FLAG_DTYPE``;
- ``units`` attributes; latitude in [-90, 90], longitude in [-180, 180];
- required dataset attributes, and that ``characterisation`` looks like
  ``<id>/<version>``;
- for every data variable: a ``<name>_flag`` variable of the same
  dimensions with the CF flag attributes of the registry, linked by
  ``ancillary_variables``; an ``origin`` attribute; and a "no value" state at
  every NaN;
- every tracked optional variable (viewing geometry, uncertainties, source
  indices ...) is either present or listed in ``availability`` with a known
  state.

Science
-------
No scientific processing.

Origin
------
New code (no MATLAB equivalent).

Validation
----------
``tests/test_data_formats.py``: a valid example of each format passes, and
each kind of fault is reported.

AI involvement
--------------
AI-originated (Claude), from the design agreed with C. Wright on 2026-11-25;
to be reviewed by C. Wright.

Usage
-----
>>> from atmosgw.curate.data import validate  # doctest: +SKIP
>>> validate(ds, "swath")  # doctest: +SKIP

Status
------
Proposal (format version 1 is a draft for agreement).

Change log
----------
2026/11/25 [Claude] initial version

How to cite
-----------
Cite the package (see ``CITATION.cff``).
"""

import re

import numpy as np

from atmosgw.curate.data import states

# variable kinds: "data" (float32, with a flag variable and an origin), "coord"
# (float64), "time" (datetime64[us]), "index" (integers or float64 indices),
# "text", "flag"
ORIGINS = ("measured", "retrieved", "derived", "mixed")

COMMON_ATTRS = (
    "format", "format_version", "instrument", "product", "characterisation",
    "data", "time_range", "availability", "source", "atmosgw_curate_version",
    "options", "files_read", "files_not_found", "multiple_matches",
    "state_counts",
)

_POINT = ("along_track", "cross_track")
_GEOMETRY = ("satellite_latitude", "satellite_longitude", "satellite_altitude",
             "view_zenith_angle", "view_azimuth_angle")

FORMATS = {
    "swath": {
        "version": 1,
        "dims": ("along_track", "cross_track", "level", "channel"),
        "attrs": COMMON_ATTRS + ("spatial_dims", "grid_spacing_km", "grid"),
        "variables": {
            "latitude": dict(kind="coord", dims=_POINT, units="degrees_north",
                             range=(-90, 90), required=True),
            "longitude": dict(kind="coord", dims=_POINT, units="degrees_east",
                              range=(-180, 180), required=True),
            "time": dict(kind="time", dims=_POINT, required=True),
            "altitude": dict(kind="coord", dims=[("level",), _POINT + ("level",)],
                             units="km", required="3d"),
            "temperature": dict(kind="data", dims=_POINT + ("level",), units="K",
                                required="3d"),
            "pressure": dict(kind="data", dims=_POINT + ("level",), units="hPa",
                             required=False),
            "brightness_temperature": dict(
                kind="data", dims=[_POINT, _POINT + ("channel",)], units="K",
                required="2d"),
            "channel_wavenumber": dict(kind="coord", dims=("channel",),
                                       units="cm-1", required=False),
            "granule": dict(kind="index", dims=("along_track",), required=True),
            "granule_day": dict(kind="time", dims=("along_track",), required=True),
            "segment": dict(kind="index", dims=("along_track",), required=True),
            "distance_along": dict(kind="coord", dims=("along_track",), units="km",
                                   required="regular"),
            "distance_cross": dict(kind="coord", dims=("cross_track",), units="km",
                                   required="regular"),
            "source_file": dict(kind="index", dims=("along_track",), tracked=True),
            "source_row": dict(kind="index", dims=("along_track",), tracked=True),
            "original_files": dict(kind="text", dims=("file",), tracked=True),
            "day_retrieval": dict(kind="index", dims=_POINT, tracked=True),
            **{name: dict(kind="coord", dims=_POINT,
                          units="km" if name == "satellite_altitude" else
                          ("degrees_north" if name == "satellite_latitude" else
                           "degrees_east" if name == "satellite_longitude" else
                           "degree"),
                          tracked=True)
               for name in _GEOMETRY},
        },
    },
    "profile": {
        "version": 1,
        "dims": ("profile", "level", "mode", "file"),
        "attrs": COMMON_ATTRS,
        "variables": {
            "latitude": dict(kind="coord", dims=("profile", "level"),
                             units="degrees_north", range=(-90, 90), required=True),
            "longitude": dict(kind="coord", dims=("profile", "level"),
                              units="degrees_east", range=(-180, 180),
                              required=True),
            "time": dict(kind="time", dims=("profile", "level"), required=True),
            "altitude": dict(kind="coord", dims=("profile", "level"), units="km",
                             required=True),
            "temperature": dict(kind="data", dims=("profile", "level"), units="K",
                                required="raw"),
            "pressure": dict(kind="data", dims=("profile", "level"), units="hPa",
                             required=True),
            "temperature_residual": dict(kind="data", dims=("profile", "level"),
                                         units="K", required="pw"),
            "temperature_pw": dict(kind="data", dims=("profile", "level", "mode"),
                                   units="K", required="pw"),
            "source_profile": dict(kind="index", dims=("profile", "level"),
                                   tracked=True),
            "source_file": dict(kind="index", dims=("profile", "level"),
                                tracked=True),
            "original_files": dict(kind="text", dims=("file",), tracked=True),
            "dx": dict(kind="coord", dims=("profile", "level"), units="km",
                       tracked=True),
        },
    },
}

_CHARACTERISATION = re.compile(r"^[A-Za-z0-9_.-]+/\d+$")


class FormatError(ValueError):
    """A Dataset does not follow its format; the message lists every problem."""


def data_variables(format_name):
    """Names of the format's data variables (those that carry flags)."""
    return [n for n, v in _spec(format_name)["variables"].items()
            if v["kind"] == "data"]


def availability(ds, name):
    """'present' if name is in ds, else its state from the 'availability' attribute.

    Raises KeyError if the variable is neither present nor listed.
    """
    if name in ds.variables:
        return "present"
    listed = states.parse_availability(ds.attrs.get("availability", ""))
    if name not in listed:
        raise KeyError(f"{name!r} is neither in the Dataset nor listed in "
                       "its 'availability' attribute.")
    return listed[name]


def validate(ds, format_name):
    """Raise FormatError listing every way ds does not follow the format."""
    found = problems(ds, format_name)
    if found:
        raise FormatError(f"Not a valid {format_name!r} Dataset:\n- "
                          + "\n- ".join(found))


def problems(ds, format_name):
    """Every way ds does not follow the format (empty list if it does)."""
    spec = _spec(format_name)
    out = []

    # identity and attributes
    if ds.attrs.get("format") != format_name:
        out.append(f"attribute 'format' is {ds.attrs.get('format')!r}, "
                   f"not {format_name!r}")
    if ds.attrs.get("format_version") != spec["version"]:
        out.append(f"attribute 'format_version' is "
                   f"{ds.attrs.get('format_version')!r}, not {spec['version']}")
    for attr in spec["attrs"]:
        if attr not in ds.attrs:
            out.append(f"missing dataset attribute {attr!r}")
    cp = ds.attrs.get("characterisation")
    if cp is not None and not _CHARACTERISATION.match(str(cp)):
        out.append(f"attribute 'characterisation' {cp!r} is not '<id>/<version>'")
    try:
        listed = states.parse_availability(ds.attrs.get("availability", ""))
    except ValueError as error:
        listed = {}
        out.append(f"attribute 'availability' cannot be read: {error}")
    for name, value in listed.items():
        if value not in states.AVAILABILITY:
            out.append(f"availability of {name!r} is {value!r}; known: "
                       f"{', '.join(states.AVAILABILITY)}")
        if name in ds.variables:
            out.append(f"{name!r} is present but also listed in 'availability'")

    # dimensions
    for dim in ds.dims:
        if dim not in spec["dims"] and dim != "file":
            out.append(f"unexpected dimension {dim!r}")

    # variables
    conditions = _conditions(ds, format_name)
    for name, var in spec["variables"].items():
        required = var.get("required", False)
        if isinstance(required, str):
            required = conditions[required]
        if name not in ds.variables:
            if required:
                out.append(f"missing required variable {name!r}")
            elif var.get("tracked") and name not in listed:
                out.append(f"optional variable {name!r} is absent but not "
                           "listed in 'availability'")
            continue
        out.extend(_check_variable(ds, name, var))
        if var["kind"] == "data":
            out.extend(_check_data_variable(ds, name, listed))
    return out


# supporting functions


def _spec(format_name):
    try:
        return FORMATS[format_name]
    except KeyError:
        raise ValueError(f"Unknown format {format_name!r}; known: "
                         f"{', '.join(FORMATS)}") from None


def _conditions(ds, format_name):
    """Which conditional requirements apply to this Dataset."""
    product = str(ds.attrs.get("product", "")).lower()
    data = str(ds.attrs.get("data", "")).lower()
    return {
        "3d": format_name == "swath" and product == "3d",
        "2d": format_name == "swath" and product == "2d",
        "regular": str(ds.attrs.get("grid", "")).lower() == "regular",
        "raw": format_name == "profile" and not data.startswith("planetary"),
        "pw": format_name == "profile" and data.startswith("planetary"),
    }


def _check_variable(ds, name, var):
    """Dimensions, type, units and range of one variable."""
    out = []
    v = ds[name]
    allowed = var["dims"] if isinstance(var["dims"], list) else [var["dims"]]
    if tuple(v.dims) not in [tuple(d) for d in allowed]:
        out.append(f"{name!r} has dimensions {tuple(v.dims)}, expected "
                   + " or ".join(str(tuple(d)) for d in allowed))
    kind, dtype = var["kind"], v.dtype
    if kind == "data" and dtype != np.float32:
        out.append(f"{name!r} is {dtype}, expected float32")
    elif kind == "coord" and dtype != np.float64:
        out.append(f"{name!r} is {dtype}, expected float64")
    elif kind == "time" and dtype != np.dtype("datetime64[us]"):
        out.append(f"{name!r} is {dtype}, expected datetime64[us]")
    if "units" in var and v.attrs.get("units") != var["units"]:
        out.append(f"{name!r} has units {v.attrs.get('units')!r}, expected "
                   f"{var['units']!r}")
    if "range" in var and v.size:
        low, high = var["range"]
        values = np.asarray(v, dtype=float)
        bad = np.isfinite(values) & ((values < low) | (values > high))
        if np.any(bad):
            out.append(f"{name!r} has {int(bad.sum())} value(s) outside "
                       f"[{low}, {high}]")
    return out


def _check_data_variable(ds, name, listed):
    """Flag variable, origin, and a 'no value' state at every NaN."""
    out = []
    v = ds[name]
    if v.attrs.get("origin") not in ORIGINS:
        out.append(f"{name!r} has origin {v.attrs.get('origin')!r}; expected one "
                   f"of {', '.join(ORIGINS)}")
    flag_name = f"{name}_flag"
    if flag_name not in str(v.attrs.get("ancillary_variables", "")).split():
        out.append(f"{name!r} does not name {flag_name!r} in ancillary_variables")
    if flag_name not in ds.variables:
        out.append(f"missing flag variable {flag_name!r}")
        return out
    f = ds[flag_name]
    if f.dims != v.dims:
        out.append(f"{flag_name!r} has dimensions {f.dims}, not those of {name!r}")
        return out
    if f.dtype != states.FLAG_DTYPE:
        out.append(f"{flag_name!r} is {f.dtype}, expected "
                   f"{np.dtype(states.FLAG_DTYPE)}")
    expected = states.flag_attrs()
    if (f.attrs.get("flag_meanings") != expected["flag_meanings"]
            or not np.array_equal(np.asarray(f.attrs.get("flag_masks", [])),
                                  expected["flag_masks"])):
        out.append(f"{flag_name!r} does not have the registry's flag_masks and "
                   "flag_meanings")
    unexplained = np.isnan(np.asarray(v, dtype=float)) & ~states.no_value_mask(f)
    if np.any(unexplained):
        out.append(f"{name!r} has {int(unexplained.sum())} NaN value(s) with no "
                   "'no value' state in its flag")
    uncertainty = f"{name}_uncertainty"
    if uncertainty not in ds.variables and uncertainty not in listed:
        out.append(f"{uncertainty!r} is absent but not listed in 'availability'")
    return out
