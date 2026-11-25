"""
Reader for the ALOMAR lidar (Andoya, Norway; yearly HDF5 files). Port of
``module_load_ALOMAR.m`` (uobGravWaves/MultiInstrumentFunctions, commit 78eef78,
with the height fix of commit e13e247, branch fix-alomar-heights).

Files: ``<ALOMAR folder>/<yyyy>_V2026-07-27_file4juli.h5``, one per year, with
one group per measurement day (groups read in name order, which is date
order). Each group holds ``T1_lidar`` (temperature, the higher-resolution
product) and ``T1_lidar_unc`` as [times x heights], ``time_lidar`` (POSIX
seconds) and ``z_lidar`` (heights, km). Each profile has its own day's heights:
the grid can change between days (in 2007, from 2007-05-14). Before the fix
(Q7 in ``docs/quirks.md``), the first day's heights were used for all days.

Profiles are the lidar's time steps. All are at the site (69.38 N, 16.08 E).
``SourceProf`` counts profiles through each yearly file in day order, and
``SourceFile`` counts the yearly files read. Pressure is estimated from height
with ``height_to_pressure``.

Days are concatenated, padded with NaN, then profiles in the (expanded) time
range are kept. Unlike the other readers, ``additional_vars`` are not passed
through ``conform_vars`` (as in MATLAB): one value per time step is repeated
over heights, anything else is used as it is.
"""

import os

import h5py
import numpy as np

from atmosgw.curate.atmosphere import height_to_pressure
from atmosgw.curate.matlab_time import datenum_to_datetime
from atmosgw.curate.obs._limb_tools import ESTIMATED_PRES, note, note_sources, record
from atmosgw.curate.profiles import concat_profiles

SITE_LAT, SITE_LON = 69.38, 16.08
RESOLUTION = "T1"  # the higher-resolution product
DATA_VERSION = "V2026-07-27"  # supersedes V2025-07-09


def load(settings, inst_info, var_names):
    """Read the years touched; returns (data, file_list). See limb.py."""
    h5_vars = list(dict.fromkeys([f"{RESOLUTION}_lidar", f"{RESOLUTION}_lidar_unc",
                                  "time_lidar", "z_lidar",
                                  *settings["additional_vars"]]))
    first, last = (int(str(datenum_to_datetime(t))[:4])
                   for t in (min(settings["time_range"]), max(settings["time_range"])))

    note_sources(inst_info, {"Temp": f"{RESOLUTION}_lidar", "Time": "time_lidar",
                             "Alt": "z_lidar"})
    note(inst_info, "Lat", "site latitude (fixed)")
    note(inst_info, "Lon", "site longitude (fixed)")
    note(inst_info, "Pres", ESTIMATED_PRES)

    days, file_list = [], []
    for year in range(first, last + 1):
        file_name = os.path.join(inst_info["path"],
                                 f"{year}_{DATA_VERSION}_file4juli.h5")
        if not os.path.isfile(file_name):
            record(inst_info)["files_not_found"].append(str(year))
            continue
        file_list.append(os.path.basename(file_name))
        profiles_in_file = 0
        with h5py.File(file_name, "r") as f:
            for group in sorted(f.keys()):
                a = {name: np.asarray(f[group][name][()], dtype=np.float64)
                     for name in h5_vars}
                if settings["verbose"]:
                    print(f"----> Loading ALOMAR data: {group}")
                # POSIX seconds -> datenum (as datenum(datetime(...,'posixtime')))
                a["time_lidar"] = a["time_lidar"] / 86400 + 719529
                n_time = a["time_lidar"].size
                a["PN"] = profiles_in_file + np.arange(1.0, n_time + 1)
                a["FN"] = np.full(n_time, float(len(file_list)))
                profiles_in_file += n_time
                # each profile's heights, from its own day (Q7, fixed)
                a["ZP"] = np.tile(a["z_lidar"][None, :], (n_time, 1))
                days.append({k: (v[:, None] if v.ndim == 1 else v)
                             for k, v in a.items() if k != "z_lidar"})

    if not days:  # no file for the years asked for: no profiles
        return {name: np.empty((0, 0)) for name in var_names}, file_list

    store = concat_profiles(days, [k for k in days[0]])
    time = store["time_lidar"][:, 0]
    keep = np.flatnonzero((time >= settings["time_range"][0])
                          & (time <= settings["time_range"][1]))
    store = {k: v[keep] for k, v in store.items()}

    temp = store[f"{RESOLUTION}_lidar"]
    n_time, n_z = temp.shape
    data = {
        "Temp": temp,
        "Lat": np.full(temp.shape, SITE_LAT),
        "Lon": np.full(temp.shape, SITE_LON),
        "Time": np.tile(store["time_lidar"], (1, n_z)),
        "Alt": store["ZP"],
    }
    data["Pres"] = height_to_pressure(data["Alt"])
    data["SourceFile"] = np.tile(store["FN"], (1, n_z))
    data["SourceProf"] = np.tile(store["PN"], (1, n_z))
    for name in settings["additional_vars"]:
        value = store[name]
        if value.shape[1] == 1 and value.shape[0] == n_time:
            value = np.tile(value, (1, n_z))  # one value per time step
        data[name] = value
    return data, file_list
