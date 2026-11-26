"""Common data formats, data states, and their checker.

``validate(ds, format_name)``: check a Dataset against a format (``"swath"``,
``"profile"``, ``"time_height"``); ``availability(ds, name)``: whether a
variable is present, or why not; ``to_profiles(ds)``: a 3-D swath or a
time_height series as profiles; ``to_netcdf(ds, path)``, ``load_netcdf(path)``:
save a valid Dataset to a netCDF file and read it back. See ``formats.py``,
``states.py``, ``convert.py``, ``netcdf.py`` and
``docs/design/data_formats.md``.
"""

from atmosgw.curate.data.convert import to_profiles
from atmosgw.curate.data.formats import (
    FORMATS,
    FormatError,
    availability,
    problems,
    validate,
)
from atmosgw.curate.data.netcdf import load_netcdf, to_netcdf

__all__ = ["FORMATS", "FormatError", "availability", "load_netcdf", "problems",
           "to_netcdf", "to_profiles", "validate"]
