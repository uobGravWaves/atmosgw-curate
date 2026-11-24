"""
Limb-sounder reader for MIPAS (Oxford retrieval, netCDF from CEDA). Port of
``module_load_MIPAS.m`` (uobGravWaves/MultiInstrumentFunctions, commit 78eef78).

Files: one per day, the first below the MIPAS folder whose path contains
``_<yyyy><mm><dd>_nom``. Variables are read as ``rCDF`` reads them (see
``readers.read_netcdf``): float64, with ``_FillValue`` (-999) as NaN.
Temperature ``TEM``, pressure ``PRE`` and height ``hgt`` are [profiles x
levels]; ``latitude``, ``longitude`` and ``time`` (hours since the start of
the file's day) are per profile.
"""

import numpy as np

from atmosgw.curate.matlab_time import datenum_to_datetime
from atmosgw.curate.obs._limb_tools import conform_vars, loop_days
from atmosgw.curate.readers import find_files, read_netcdf

# output names -> MIPAS variable names
FILE_NAMES = {"Temp": "TEM", "Lat": "latitude", "Lon": "longitude", "Alt": "hgt",
              "Pres": "PRE", "Time": "time"}


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""

    def find_file(day):
        date = str(datenum_to_datetime(day))[:10].replace("-", "")
        found = find_files(inst_info["path"], f"_{date}_nom")
        return found[0] if found else None

    def read_file(path, file_count, day):
        names = [FILE_NAMES.get(n, n) for n in var_names
                 if n not in ("SourceProf", "SourceFile")]
        try:
            working = read_netcdf(path, names)
        except KeyError as error:
            raise KeyError(f"Variable {error} not found in {path}") from error
        store = {}
        for name in var_names:
            if name == "Time":
                store["Time"] = day + _column(working["time"]) / 24
            elif name == "SourceProf":
                n = store["Lat"].size
                store["SourceProf"] = np.arange(1.0, n + 1)[:, None]  # 1-based
            elif name == "SourceFile":
                store["SourceFile"] = np.ones_like(store["SourceProf"]) * file_count
            else:
                store[name] = _column(working[FILE_NAMES.get(name, name)])

        n_lev = store["Alt"].shape[1]
        for name in ("Lat", "Lon", "Time", "SourceProf", "SourceFile"):
            store[name] = np.tile(store[name], (1, n_lev))
        return conform_vars(store, settings)

    return loop_days(settings, var_names, find_file, read_file, "MIPAS")


def _column(value):
    """1-D variables as columns, as rCDF gives them."""
    return value[:, None] if value.ndim == 1 else value
