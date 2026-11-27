"""
Limb-sounder reader for SOFIE (GATS Level 2 netCDF files). Port of
``module_load_SOFIE.m`` (uobGravWaves/MultiInstrumentFunctions, commit 78eef78).

Files: one per day, the file below the SOFIE folder whose path contains
``Level2_<yyyy><ddd>`` (of several, the newest version in the name, e.g.
``_01.3.nc``, else the first: Q19). Variables are read as ``rCDF`` reads them (see
``readers.read_netcdf``): float64, with ``missing_value`` (-1e24) as NaN.
``Temperature`` and ``Pressure`` are [events x altitudes]; ``Altitude`` is one
set of levels for all events; positions are those at 83 km
(``Latitude_83km``, ``Longitude_83km``). ``Time_UT`` is in hours since the
start of the day, although the files give its units as "msec" (checked on
2007-06-01: values 0.87-23.4).
"""

import numpy as np

from atmosgw.curate.obs._common import (
    choose_file,
    day_label,
    find_in,
    note,
    note_sources,
    note_version,
    year_and_day,
)
from atmosgw.curate.obs._limb_tools import (
    VERSION_IN_NAME,
    conform_vars,
    loop_days,
    newest_version,
)
from atmosgw.curate.readers import read_attrs, read_netcdf

ENDORSED_KM = (15.0, 88.0)  # temperature (SOFIE V1.3 data description, Table 1)

# output names -> SOFIE variable names
FILE_NAMES = {"Temp": "Temperature", "Lat": "Latitude_83km", "Lon": "Longitude_83km",
              "Alt": "Altitude", "Pres": "Pressure", "Time": "Time_UT"}


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    note_sources(inst_info, dict(FILE_NAMES))
    note(inst_info, "Time", "file Time_UT: hours since the start of the day "
                            "(the file's units attribute says msec)")
    note(inst_info, "Lat", "position at 83 km, used for the whole profile")
    note(inst_info, "Lon", "position at 83 km, used for the whole profile")

    def find_file(day):
        year, doy = year_and_day(day)
        found = find_in(inst_info, f"Level2_{year:04d}{doy:03d}")
        # several versions of a day: the newest (Q19; MATLAB takes the first)
        return choose_file(found, inst_info, day_label(day),
                           newest_version(VERSION_IN_NAME["SOFIE"]))

    def read_file(path, file_count, day):
        names = [FILE_NAMES.get(n, n) for n in var_names
                 if n not in ("SourceProf", "SourceFile", "InstrumentBad")]
        if "InstrumentBad" in var_names and "Altitude" not in names:
            names.append("Altitude")
        try:
            working = read_netcdf(path, names)
        except KeyError as error:
            raise KeyError(f"Variable {error} not found in {path}") from error
        version = read_attrs(path).get("DP_Version", "")
        note_version(inst_info, path, f"DP_Version {version}" if version else "")
        store = {}
        for name in var_names:
            if name == "InstrumentBad":
                # heights not endorsed for scientific use (SOFIE V1.3 data
                # description, Table 1: temperature endorsed 15-88 km, reported
                # 10-102 km; checked 2026-11-26)
                alt = working["Altitude"][None, :]
                with np.errstate(invalid="ignore"):
                    bad = (alt < ENDORSED_KM[0]) | (alt > ENDORSED_KM[1])
                store[name] = np.tile(bad, (working["Temperature"].shape[0]
                                            if "Temperature" in working else 1, 1)
                                      ).astype(float)
                note(inst_info, "Temp", "heights outside 15-88 km (not endorsed for "
                                        "scientific use, SOFIE V1.3 data "
                                        "description) are instrument_flag")
            elif name == "Time":
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

    return loop_days(settings, var_names, find_file, read_file, "SOFIE", inst_info)
