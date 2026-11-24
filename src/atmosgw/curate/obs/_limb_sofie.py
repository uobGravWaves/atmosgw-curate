"""
Limb-sounder reader for SOFIE (GATS Level 2 netCDF files). Port of
``module_load_SOFIE.m`` (uobGravWaves/MultiInstrumentFunctions, commit 78eef78).

Files: one per day, the first below the SOFIE folder whose path contains
``Level2_<yyyy><ddd>``. Variables are read as ``rCDF`` reads them (see
``readers.read_netcdf``): float64, with ``missing_value`` (-1e24) as NaN.
``Temperature`` and ``Pressure`` are [events x altitudes]; ``Altitude`` is one
set of levels for all events; positions are those at 83 km
(``Latitude_83km``, ``Longitude_83km``). ``Time_UT`` is in hours since the
start of the day, although the files give its units as "msec" (checked on
2007-06-01: values 0.87-23.4).
"""

import numpy as np

from atmosgw.curate.obs._limb_tools import conform_vars, loop_days, year_and_day
from atmosgw.curate.readers import find_files, read_netcdf

# output names -> SOFIE variable names
FILE_NAMES = {"Temp": "Temperature", "Lat": "Latitude_83km", "Lon": "Longitude_83km",
              "Alt": "Altitude", "Pres": "Pressure", "Time": "Time_UT"}


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""

    def find_file(day):
        year, doy = year_and_day(day)
        found = find_files(inst_info["path"], f"Level2_{year:04d}{doy:03d}")
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
                store["Time"] = day + working["Time_UT"][:, None] / 24
            elif name == "Alt":
                store["Alt"] = working["Altitude"][None, :]  # one row of levels
            elif name == "SourceProf":
                n_prof, n_lev = store["Lat"].shape[0], store["Alt"].shape[1]
                profile = np.arange(1.0, n_prof + 1)[:, None]  # 1-based
                store["SourceProf"] = np.tile(profile, (1, n_lev))
            elif name == "SourceFile":
                store["SourceFile"] = np.ones_like(store["SourceProf"]) * file_count
            else:
                value = working[FILE_NAMES.get(name, name)]
                store[name] = value[:, None] if value.ndim == 1 else value

        n_lev, n_prof = store["Alt"].size, store["Lat"].size
        for name in ("Lat", "Lon", "Time"):
            store[name] = np.tile(store[name], (1, n_lev))
        store["Alt"] = np.tile(store["Alt"], (n_prof, 1))
        return conform_vars(store, settings)

    return loop_days(settings, var_names, find_file, read_file, "SOFIE")
