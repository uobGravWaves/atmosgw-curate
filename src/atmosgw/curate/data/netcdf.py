"""
Write Datasets in the common formats to netCDF files, and read them back.

``to_netcdf`` saves a loader's output (any of the formats in
``data.formats``) so that ``load_netcdf`` returns it as it was: the same
variables, types, states and attributes, with times exact to the
microsecond. Only valid Datasets are written, and a failed write leaves no
file behind, so a file that exists was written completely from a valid
Dataset.

Inputs and outputs
------------------
``to_netcdf(ds, path, pack=False, clobber=False)`` -> None

``load_netcdf(path)`` -> ``xarray.Dataset``

===========  ===========  =======  ==============================================
Option       Type         Default  Description
===========  ===========  =======  ==============================================
``ds``       Dataset      -        a Dataset in one of the formats ("profile",
                                   "swath", "time_height", "trajectory"), as the
                                   loaders return
``path``     str or Path  -        the file to write (its folder must exist)
``pack``     bool         False    store suitable variables as scaled integers
                                   (lossy, to a stated resolution; see below)
``clobber``  bool         False    replace an existing file (otherwise an
                                   existing file is an error)
===========  ===========  =======  ==============================================

``to_netcdf``:

- checks ``ds`` with ``data.validate`` against the format its ``format``
  attribute names, and raises ``FormatError`` (writing nothing) if it is not
  valid;
- writes netCDF-4 with zlib compression (level 4, shuffle) of every numeric
  variable; nothing else changes the values unless ``pack``;
- stores times (``time``, ``granule_day``, ``native_time``) as int64
  microseconds since 1970-01-01, with NaT as ``_FillValue``, so they read back
  as the same ``datetime64[us]`` values;
- writes to a temporary file in the same folder and renames it to ``path``
  only when complete; on any error the temporary file is removed and ``path``
  is left as it was;
- ignores any encoding the Dataset carries from a file it was read from (so
  re-saving a packed file without ``pack`` stores the values unpacked);
- text variables cannot hold NaN in netCDF: missing values in text (e.g.
  some GNSS ``MetaData`` columns loaded with ``additional_vars``) are stored,
  and read back, as empty text.

``pack=True`` stores these variables as integers with ``scale_factor`` and
``add_offset`` (CF packing) at the resolution given, ``_FillValue`` for NaN:

=====================================================  ===============  ==========
Variables                                              Resolution       Why
=====================================================  ===============  ==========
data variables in K (``temperature``,                  0.01 K           10x finer
``temperature_residual``, ``temperature_pw``,                           than the
``brightness_temperature``,                                             best data
``brightness_temperature_perturbation``)                                (0.1 K)
their uncertainties in K (``<name>_uncertainty``)      0.001 K
``latitude``, ``longitude``, ``native_latitude``,      1e-7 degree      about 1 cm
``native_longitude``
``altitude``, ``dx``                                   1e-6 km (1 mm)
=====================================================  ===============  ==========

Positions and heights whose values are all exactly float32 values (read
from float32 files, not computed: e.g. native SABER or AIRS positions) are
left unpacked: they compress well as they are, and packing made such files
larger (native SABER 2015-07-01: 5.3 MB with temperature packed, 6.1 MB with
positions packed too). Each packed variable is stored as int16 if its values
fit at that resolution, else int32, else unpacked; a value is never
clipped, so every value read back is within half the resolution of the value
written (plus the rounding of the variable's own floating-point type).
Variables with fewer than two dimensions (grid axes such as a swath's
``altitude`` or ``distance_along``, the site of a
time_height series), variables with no finite value or with infinities,
pressure (it spans several orders of magnitude), times, flags, indices, text
and ``additional_vars`` are never packed. Packed variables read back in their
own type (float32 or float64), as the formats require. Packing is lossy: use
it for archives, not for files that must reproduce a result bit for bit.

``load_netcdf``: reads a file written by ``to_netcdf`` (or any netCDF file
of a Dataset in one of the formats) into memory, unpacks it, decodes times
as ``datetime64[us]`` and checks the result with ``data.validate``
(``FormatError`` if it is not valid; to read such a file anyway, use
``xarray.load_dataset``).

Origin
------
New code (no MATLAB equivalent), from the design in
``docs/design/data_formats.md`` section 11 (agreed with C. Wright,
2026-11-25), with packing of positions and heights and the names ``clobber``
and ``load_netcdf`` agreed with C. Wright on 2026-11-26.

Validation
----------
``tests/test_netcdf.py``: each format written and read back unpacked is
identical (values, types, NaN and NaT, flags, attributes); packed values
within half their resolution and of their own type; integer type chosen by
range, unpacked when out of range; float32-exact positions left unpacked;
invalid Datasets, existing files and failed writes leave no file; re-saving
a packed file unpacked; empty Datasets. Real data (2026-11-26, eepc-0184):
SABER 2007-01-02 and 2015-07-01 (native), GNSS 2007-01-02 with MetaData,
AIRS and IASI two and one hours, ALOMAR 2007-01-23: unpacked identical
(except missing MetaData text, read back empty), packed within half the
resolution; files 2-4 times smaller than xarray's default unpacked, about
half again packed.

Status
------
Proposal. AI-originated (Claude), from the design agreed with C. Wright; to
be reviewed by C. Wright.

Usage
-----
>>> from atmosgw.curate.data import load_netcdf, to_netcdf  # doctest: +SKIP
>>> to_netcdf(ds, "/path/to/saber_2007-01-01.nc")  # doctest: +SKIP
>>> back = load_netcdf("/path/to/saber_2007-01-01.nc")  # doctest: +SKIP

Change log
----------
2026/11/26 [Claude] initial version

"""

import os
import uuid
from pathlib import Path

import numpy as np
import xarray as xr

from atmosgw.curate.data.formats import FORMATS, FormatError, data_variables, validate

# resolution of packed values (pack=True), by variable kind
TEMPERATURE_RESOLUTION = 0.01  # K, data variables in K
UNCERTAINTY_RESOLUTION = 0.001  # K, their uncertainties
POSITION_RESOLUTION = 1e-7  # degrees (about 1 cm)
HEIGHT_RESOLUTION = 1e-6  # km (1 mm)

_POSITIONS = ("latitude", "longitude", "native_latitude", "native_longitude")
_HEIGHTS = ("altitude", "dx")
_TIME_UNITS = "microseconds since 1970-01-01"
_COMPRESSION = {"zlib": True, "complevel": 4, "shuffle": True}


def to_netcdf(ds, path, pack=False, clobber=False):
    """Write a Dataset in one of the formats to netCDF. See the module docstring."""
    if not isinstance(pack, (bool, np.bool_)) or not isinstance(clobber,
                                                                (bool, np.bool_)):
        raise TypeError("pack and clobber must be True or False.")
    format_name = ds.attrs.get("format")
    if format_name not in FORMATS:
        raise FormatError(
            f"The Dataset's 'format' attribute is {format_name!r}, not one of "
            f"{', '.join(FORMATS)}: to_netcdf writes only Datasets in the "
            "package's formats (use the Dataset's own to_netcdf for others).")
    validate(ds, format_name)  # (raises FormatError listing every problem)

    path = Path(os.path.expandvars(str(path))).expanduser()
    if not path.parent.is_dir():
        raise FileNotFoundError(f"The folder {path.parent} does not exist.")
    if path.exists() and not clobber:
        raise FileExistsError(f"{path} exists: pass clobber=True to replace it.")

    # (a shallow copy without the encodings of any file it was read from)
    out = ds.copy(deep=False)
    for variable in out.variables.values():
        variable.encoding = {}
    encoding = {name: _encoding(name, out[name], format_name, pack)
                for name in out.variables}

    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex[:12]}.partial")
    try:
        out.to_netcdf(temporary, engine="netcdf4", format="NETCDF4",
                      encoding=encoding)
        if path.exists() and not clobber:  # (created by someone else meanwhile)
            raise FileExistsError(f"{path} exists: pass clobber=True to replace it.")
        os.replace(temporary, path)
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise


def load_netcdf(path):
    """Read a file written by to_netcdf, checked. See the module docstring."""
    path = Path(os.path.expandvars(str(path))).expanduser()
    ds = xr.load_dataset(path, engine="netcdf4",
                         decode_times=xr.coders.CFDatetimeCoder(time_unit="us"))
    format_name = ds.attrs.get("format")
    if format_name not in FORMATS:
        raise FormatError(f"{path}: the 'format' attribute is {format_name!r}, not "
                          f"one of {', '.join(FORMATS)} (read it with "
                          "xarray.load_dataset instead).")
    validate(ds, format_name)
    return ds


# supporting functions


def _encoding(name, variable, format_name, pack):
    """netCDF encoding of one variable: compression, time units, packing."""
    kind = variable.dtype.kind
    if kind == "M":  # times: microseconds, exact
        return {"units": _TIME_UNITS, "dtype": "int64",
                "_FillValue": np.iinfo(np.int64).min, **_compressed(variable)}
    if kind not in "biuf":  # (text)
        return {}
    encoding = dict(_compressed(variable))
    resolution = _resolution(name, variable, format_name) if pack else None
    if resolution is not None:
        encoding.update(_packing(variable, resolution))
    return encoding


def _compressed(variable):
    """Compression settings, none for empty variables (no chunks to compress)."""
    return dict(_COMPRESSION) if variable.size else {}


def _resolution(name, variable, format_name):
    """The packing resolution of a variable, or None if it is never packed."""
    if variable.ndim < 2 or variable.dtype.kind != "f":
        return None
    data_names = data_variables(format_name)
    in_kelvin = variable.attrs.get("units") == "K"
    if name in data_names and in_kelvin:
        return TEMPERATURE_RESOLUTION
    if (name.endswith("_uncertainty") and name[: -len("_uncertainty")] in data_names
            and in_kelvin):
        return UNCERTAINTY_RESOLUTION
    if name in _POSITIONS + _HEIGHTS:
        # values read from float32 files (not computed) already compress
        # well, and packing them would only make the file larger and lose
        # their last digits (native SABER 2015-07-01: 5.3 MB, 6.1 MB packed)
        if _float32_exact(variable):
            return None
        return POSITION_RESOLUTION if name in _POSITIONS else HEIGHT_RESOLUTION
    return None


def _float32_exact(variable):
    """True if every value is exactly a float32 value."""
    values = np.asarray(variable.values)
    if values.dtype == np.float32:
        return True
    with np.errstate(over="ignore", invalid="ignore"):
        return bool(np.array_equal(values.astype(np.float32).astype(values.dtype),
                                   values, equal_nan=True))


def _packing(variable, resolution):
    """CF packing of a variable at resolution: the smallest of int16 and
    int32 that holds every value without clipping ({} if neither does, or
    there is nothing to pack). Scale and offset have the variable's own
    floating-point type, so it reads back in that type."""
    values = np.asarray(variable.values, dtype=np.float64)
    if np.any(np.isinf(values)):
        return {}
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {}
    low, high = float(finite.min()), float(finite.max())
    offset = resolution * np.round((low + high) / 2 / resolution)
    steps = np.round((np.array([low, high]) - offset) / resolution)
    real = variable.dtype.type
    for int_type in (np.int16, np.int32):
        info = np.iinfo(int_type)
        if info.min < steps[0] and steps[1] <= info.max:  # (the minimum is the fill)
            return {"dtype": np.dtype(int_type).name, "scale_factor": real(resolution),
                    "add_offset": real(offset), "_FillValue": int_type(info.min)}
    return {}
