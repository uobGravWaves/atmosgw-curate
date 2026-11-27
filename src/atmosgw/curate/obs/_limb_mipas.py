"""
Limb-sounder reader for MIPAS (Oxford retrieval, netCDF from CEDA). Port of
``module_load_MIPAS.m`` (uobGravWaves/MultiInstrumentFunctions, commit 78eef78).

Files: one per day, the file below the MIPAS folder whose path contains
``_<yyyy><mm><dd>_nom`` (of several, the newest version in the name, e.g.
``_v1.40.nc``, else the first: Q19). Variables are read as ``rCDF`` reads
them (see ``readers.read_netcdf``): float64, with ``_FillValue`` (-999) as
NaN.
Temperature ``TEM``, pressure ``PRE`` and height ``hgt`` are [profiles x
levels]; ``latitude``, ``longitude`` and ``time`` (hours since the start of
the file's day) are per profile.
"""

import os
import re

import numpy as np

from atmosgw.curate.matlab_time import datenum_to_datetime
from atmosgw.curate.obs._common import (
    choose_file,
    day_label,
    find_in,
    note,
    note_sources,
    note_version,
)
from atmosgw.curate.obs._limb_tools import (
    VERSION_IN_NAME,
    conform_vars,
    loop_days,
    newest_version,
)
from atmosgw.curate.readers import read_attrs, read_netcdf

# output names -> MIPAS variable names
FILE_NAMES = {"Temp": "TEM", "Lat": "latitude", "Lon": "longitude", "Alt": "hgt",
              "Pres": "PRE", "Time": "time"}


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    note_sources(inst_info, dict(FILE_NAMES))
    note(inst_info, "Time", "file time: hours since the start of the file's day")

    def find_file(day):
        date = str(datenum_to_datetime(day))[:10].replace("-", "")
        found = find_in(inst_info, f"_{date}_nom")
        # several versions of a day: the newest (Q19; MATLAB takes the first)
        return choose_file(found, inst_info, day_label(day),
                           newest_version(VERSION_IN_NAME["MIPAS"]))

    def read_file(path, file_count, day):
        names = [FILE_NAMES.get(n, n) for n in var_names
                 if n not in ("SourceProf", "SourceFile")]
        try:
            working = read_netcdf(path, names)
        except KeyError as error:
            raise KeyError(f"Variable {error} not found in {path}") from error
        note_version(inst_info, path, _version(path))
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
            if name not in store:  # (source indices only if asked for)
                continue
            store[name] = np.tile(store[name], (1, n_lev))
        return conform_vars(store, settings)

    return loop_days(settings, var_names, find_file, read_file, "MIPAS", inst_info)


def _version(path):
    """The L2 processing version in the file's history, else in its name."""
    found = re.search(r"L2 processing version ([0-9.]+)",
                      read_attrs(path).get("history", ""))
    if found:
        return f"L2 processing version {found.group(1)}"
    found = re.search(r"_v([0-9.]+)\.nc$", os.path.basename(path))
    return f"v{found.group(1)} (file name)" if found else ""


def _column(value):
    """1-D variables as columns, as rCDF gives them."""
    return value[:, None] if value.ndim == 1 else value
