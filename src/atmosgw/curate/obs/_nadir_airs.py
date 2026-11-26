"""
Nadir-sounder reader for AIRS 3-D temperature retrievals (Hoffmann and
Alexander), in the group's netCDF files.

Port of the AIRS file handling of ``prep_airs_3d.m`` (its
``get_airs_granule``; corwin365/MatlabFunctions, commit f1a78c0).

Files, below the AIRS folder (``datasets.toml``, ``AIRS`` ``raw``), as MATLAB:

- daily files ``YYYY/airs_YYYYdDDD.nc``, holding all 240 granules of a day
  (dimension ``time``), each 135 along-track x 90 cross-track points on 27
  levels; used if present;
- otherwise granule files ``YYYY/DDD/airs_YYYY_DDD_GGG.nc``.

Granules: 240 a day, each six minutes; granule g of a day nominally starts at
00:05:28 + 6(g - 1) minutes UTC (2007 files; used only to find candidates,
with a 15-minute margin). Each granule's actual time is the mean of its
``l1_time``; a granule belongs to the time range if that mean lies in it (as
``load_limb_sounders`` selects whole profiles by their mean time).

Variables (reader contract, ``docs/design/nadir_loader.md`` section 5):
``Lat``, ``Lon``, ``Time`` (datenums), ``Temp`` [cross x along (x level)],
``Alt`` [level] (km), and per along-track row ``Granule``, ``Day``
(datenum of the granule's day), ``SourceFile`` (1-based, into the file list)
and ``SourceRow`` (1-based row within its granule). Time is converted from
seconds since 2000-01-01T00:00Z, leap seconds ignored (as MATLAB; agreed).
If a granule's heights are all zero or NaN, the standard 27 AIRS levels are
used, as MATLAB (``ret_z`` missing in some granules), and noted.

Fill values (``_common.nadir_fill_to_nan``; the files declare none): -999 (in
padded daily-file granules), -9999 positions (missing scans, times kept) and
0/0 positions or zero times become NaN; stray pixel times (> 600 s from their
scan) and scans more than a day outside the file's day (e.g. 2009 days 193-194)
become NaN. Counts are in the load record's notes. A granule left with no
valid time (e.g. all -999) is skipped and recorded in ``files_not_found`` as
"(no valid times)"; for a time range it is placed by its nominal time.

Differences from MATLAB, by design: a missing granule is skipped and recorded
(MATLAB fails the whole call); granules are returned in time order; the fill
handling above (MATLAB reads fill values as data, and crashes regridding them).

Note: the daily files on eepc-0184 store temperature rounded to whole kelvin
(``merge_airs_into_daily.sh``); the granule files have full precision. Both
are read the same way (``readers.read_netcdf`` applies any scale factor).

Change log:
2026/11/26 [Claude] fill values, bad times to NaN; fill granules recorded, not dropped
"""

import os

import numpy as np

from atmosgw.curate.obs._common import (
    day_label,
    nadir_fill_to_nan,
    note,
    note_fill_counts,
    note_sources,
    record,
    year_and_day,
)
from atmosgw.curate.readers import read_netcdf

N_CROSS = 90
ROWS_PER_GRANULE = 135
GRANULES_PER_DAY = 240
GRANULE_DAYS = 6 / 1440  # six minutes
FIRST_START = (5 * 60 + 28) / 86400  # nominal start of granule 1, as a day fraction
J2000 = 730486.0  # MATLAB datenum(2000,1,1)

# the 27 standard AIRS retrieval levels (km), used when a granule has none (MATLAB)
AIRS_LEVELS = np.array([0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45,
                        48, 51, 54, 57, 60, 65, 70, 75, 80, 85, 90], dtype=float)

FILE_NAMES = {"Lat": "l1_lat", "Lon": "l1_lon", "Temp": "ret_temp", "Alt": "ret_z"}


def candidate_granules(settings):
    """(day, granule) pairs that may hold data in the time range, or those asked for."""
    if settings.get("granules") is not None:
        return list(settings["granules"])
    t0, t1 = settings["time_range"]
    margin = 15 / 1440
    out = []
    for day in range(int(np.floor(t0)) - 1, int(np.floor(t1)) + 1):
        for g in range(1, GRANULES_PER_DAY + 1):
            mid = day + FIRST_START + (g - 0.5) * GRANULE_DAYS
            if t0 - margin <= mid <= t1 + margin:
                out.append((day, g))
    return out


# what was done to the daily files before loading (output attribute preprocessing)
DAILY_PREPROCESSING = ("daily files merged from L. Hoffmann's granule files "
                       "(merge_airs_into_daily.sh): temperatures rounded to whole "
                       "kelvin, times stored as float32")


def find_file(root, day, granule):
    """(path, kind): the daily file if present, else the granule file; (None, None)."""
    year, doy = year_and_day(day)
    daily = os.path.join(root, f"{year:04d}", f"airs_{year:04d}d{doy:03d}.nc")
    if os.path.isfile(daily):
        return daily, "daily"
    single = os.path.join(root, f"{year:04d}", f"{doy:03d}",
                          f"airs_{year:04d}_{doy:03d}_{granule:03d}.nc")
    if os.path.isfile(single):
        return single, "granule"
    return None, None


def load(settings, inst_info, var_names):
    """Read the granules; returns (data, file_list). See the module docstring."""
    note_sources(inst_info, {**FILE_NAMES, "Time": "l1_time"})
    note(inst_info, "Time", "l1_time: seconds since 2000-01-01T00:00Z")
    extra = [n for n in var_names if n not in ("Lat", "Lon", "Time", "Temp", "Alt",
                                                "Granule", "Day", "SourceFile",
                                                "SourceRow")]
    file_names = ["l1_time", "l1_lat", "l1_lon", "ret_temp"] + extra
    explicit = settings.get("granules") is not None
    t0, t1 = (-np.inf, np.inf) if explicit else settings["time_range"]

    # group the candidates by file, so each daily file is opened once
    by_file, not_found = {}, []
    for day, g in candidate_granules(settings):
        path, kind = find_file(inst_info["path"], day, g)
        if path is None:
            not_found.append(f"{day_label(day)} granule {g}")
            continue
        by_file.setdefault((path, kind), []).append((day, g))

    granules, file_list, alt = [], [], None
    fill_counts = {}
    for (path, kind), wanted in sorted(by_file.items()):
        rows = [g - 1 for _, g in wanted] if kind == "daily" else None
        if kind == "daily":  # (granule files are the producer's, unaltered)
            done = record(inst_info).setdefault("preprocessing", set())
            done.add(DAILY_PREPROCESSING)
        values = read_netcdf(path, file_names, rows=rows)
        z = read_netcdf(path, ["ret_z"])["ret_z"].ravel()
        if np.nansum(np.abs(z)) == 0:
            z = AIRS_LEVELS
            note(inst_info, "Alt", "file heights missing: standard 27 AIRS levels used")
        used = False
        for k, (day, g) in enumerate(wanted):
            part = {name: (v[k] if kind == "daily" else v)
                    for name, v in values.items()}
            span = ((day - J2000) * 86400, (day + 1 - J2000) * 86400)
            counts = nadir_fill_to_nan(part["l1_lat"], part["l1_lon"],
                                       part["l1_time"], part["ret_temp"], span)
            time = J2000 + part["l1_time"].T / 86400  # [cross x along], datenums
            mean = np.nanmean(time) if np.any(np.isfinite(time)) else np.nan
            if not np.isfinite(mean):  # in the range by its nominal time?
                nominal = day + FIRST_START + (g - 0.5) * GRANULE_DAYS
                if explicit or t0 <= nominal <= t1:
                    not_found.append(f"{day_label(day)} granule {g} (no valid times)")
                continue
            if not explicit and not (t0 <= mean <= t1):
                continue
            for key, n in counts.items():  # (granules returned only)
                fill_counts[key] = fill_counts.get(key, 0) + n
            if not used:
                file_list.append(os.path.basename(path))
                used = True
            if alt is None:
                alt = z
            elif not np.array_equal(alt, z):
                raise ValueError(f"{os.path.basename(path)}: heights differ from "
                                 "those of the granules before it.")
            n_along = part["l1_lat"].shape[0]
            g_data = {
                "Lat": part["l1_lat"].T, "Lon": part["l1_lon"].T, "Time": time,
                "Temp": np.transpose(part["ret_temp"], (1, 0, 2)),
                "Granule": np.full(n_along, float(g)),
                "Day": np.full(n_along, float(day)),
                "SourceFile": np.full(n_along, float(len(file_list))),
                "SourceRow": np.arange(1, n_along + 1, dtype=float),
            }
            for name in extra:
                v = np.asarray(part[name])
                g_data[name] = (np.transpose(v, (1, 0) + tuple(range(2, v.ndim)))
                                if v.ndim >= 2 else v)
            granules.append((mean, g_data))
    record(inst_info)["files_not_found"].extend(not_found)
    note_fill_counts(inst_info, fill_counts)

    granules.sort(key=lambda item: item[0])
    names = ["Lat", "Lon", "Time", "Temp", "Granule", "Day", "SourceFile",
             "SourceRow"] + extra
    if not granules:
        n_lev = AIRS_LEVELS.size
        data = {"Lat": np.empty((N_CROSS, 0)), "Lon": np.empty((N_CROSS, 0)),
                "Time": np.empty((N_CROSS, 0)), "Temp": np.empty((N_CROSS, 0, n_lev)),
                **{n: np.empty(0) for n in ("Granule", "Day", "SourceFile",
                                            "SourceRow")}}
        data["Alt"] = AIRS_LEVELS
        return data, file_list
    data = {}
    for name in names:
        parts = [g[name] for _, g in granules]
        data[name] = np.concatenate(parts, axis=1 if parts[0].ndim >= 2 else 0)
    data["Alt"] = alt
    return data, file_list


def nominal_granule(t):
    """The AIRS granule (day, number) whose nominal span contains datenum t."""
    day = int(np.floor(t - FIRST_START))
    g = int(np.floor((t - day - FIRST_START) / GRANULE_DAYS)) + 1
    if g > GRANULES_PER_DAY:
        day, g = day + 1, g - GRANULES_PER_DAY
    return day, g


__all__ = ["load", "candidate_granules", "find_file", "nominal_granule",
           "AIRS_LEVELS", "N_CROSS", "ROWS_PER_GRANULE"]
