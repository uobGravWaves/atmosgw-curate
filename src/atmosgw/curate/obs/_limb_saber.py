"""
Limb-sounder reader for SABER (GATS monthly Temp_O3_H2O netCDF files, v2.0).
Port of ``module_load_SABER.m`` (uobGravWaves/MultiInstrumentFunctions,
commit 78eef78).

Files: one per month, the first below the SABER folder whose path contains
``<Month><yyyy>`` (e.g. ``January2007``, any case). Variables are read as
``rCDF`` reads them (see ``readers.read_netcdf``): float64, with
``missing_value`` (-999) as NaN. Each event (profile) has a date ``date``
(yyyyddd) and, at each level, ``time`` in ms since midnight of that date
(it can exceed a day for orbits crossing midnight).

Days: from the day before the start (files are orbit-based, so a day's
profiles can be in the previous day's month file) to the end day. Each day
takes the events whose date is that day. Altitudes above 1000 km are set to
NaN, as are negative times.

Only the variables needed are read, and for each day only that day's events;
the result is the same as MATLAB's, which reads whole files.
"""

import math
import os

import numpy as np

from atmosgw.curate.matlab_time import datenum, datenum_to_datetime
from atmosgw.curate.obs._common import (
    choose_file,
    note,
    note_sources,
    note_version,
    record,
)
from atmosgw.curate.obs._limb_tools import conform_vars
from atmosgw.curate.profiles import concat_profiles
from atmosgw.curate.readers import find_files, read_attrs, read_netcdf

MONTHS = ("January", "February", "March", "April", "May", "June", "July",
          "August", "September", "October", "November", "December")

# output names -> SABER variable names
FILE_NAMES = {"Lat": "tplatitude", "Lon": "tplongitude", "Pres": "pressure",
              "Alt": "tpaltitude", "Temp": "ktemp"}


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    first = math.floor(min(settings["time_range"])) - 1  # and the day before
    last = math.floor(max(settings["time_range"]))
    note_sources(inst_info, {**FILE_NAMES, "Time": "date, time"})
    note(inst_info, "Time", "file date (yyyyddd) + time (ms since midnight)")
    file_name, event_dates, chosen = "", None, None
    months_seen = set()  # (each month is looked up for every day, recorded once)
    days, file_list = [], []
    for k, day in enumerate(range(first, last + 1), start=1):
        if settings["verbose"]:
            print(f"----> Loading SABER data: day {k} of {last - first + 1}")
        t = datenum_to_datetime(day)
        year = int(t.astype("datetime64[Y]").astype(int)) + 1970
        month = MONTHS[int(t.astype("datetime64[M]").astype(int)) % 12]
        found = find_files(inst_info["path"], f"{month}{year}")
        first_look = f"{month}{year}" not in months_seen
        months_seen.add(f"{month}{year}")
        if first_look:  # (the file of the month, recording duplicates)
            chosen = choose_file(found, inst_info, f"{month}{year}")
        if chosen is None:
            # (the day before the request is only looked at: not reported)
            if first_look and day > first:
                record(inst_info)["files_not_found"].append(f"{month}{year}")
            continue

        if chosen != file_name:  # monthly files: read only when it changes
            file_name = chosen
            file_list.append(os.path.basename(file_name))
            version = read_attrs(file_name).get("Data_Product_Version", "")
            note_version(inst_info, file_name,
                         f"Data_Product_Version {version}" if version else "")
            date = read_netcdf(file_name, ["date"])["date"]
            year_of = np.floor(date / 1000)
            event_dates = datenum(year_of, 1, date - year_of * 1000)

        events = np.flatnonzero(np.floor(event_dates) == day)
        if events.size == 0:
            continue

        names = ["date", "time"] + [FILE_NAMES.get(n, n) for n in var_names
                                    if n not in ("Time", "SourceProf", "SourceFile")]
        working = read_netcdf(file_name, list(dict.fromkeys(names)), rows=events)
        n_lev = working["time"].shape[1]

        store = {}
        for name in var_names:
            if name == "SourceProf":
                store[name] = np.tile((events + 1.0)[:, None], (1, n_lev))  # 1-based
            elif name == "SourceFile":
                store[name] = np.full((events.size, n_lev), float(len(file_list)))
            elif name != "Time":
                store[name] = working[FILE_NAMES.get(name, name)]

        # bad heights break the interpolation
        store["Alt"] = np.where(store["Alt"] > 1000, np.nan, store["Alt"])

        # time: ms since midnight of each event's date (may exceed one day)
        date = np.tile(working["date"][:, None], (1, n_lev))
        year_of = np.floor(date / 1000)
        day_of_year = date - year_of * 1000
        seconds = working["time"] / 1000
        seconds[seconds < 0] = np.nan
        store["Time"] = datenum(year_of, 1, day_of_year, 0, 0, seconds)

        days.append(conform_vars(store, settings))

    return concat_profiles(days, var_names, consume=True), file_list
