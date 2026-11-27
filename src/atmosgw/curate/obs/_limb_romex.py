"""
Limb-loader reader for ROMEX GNSS radio-occultation profiles, 2022-09-01 to
2022-11-30 (ROMEX taken to be the IROWG Radio Occultation Modeling
Experiment: UNVERIFIED).

Files: one or more per mission, ``romex_<mission>rmx[...]_100m_v1a[...].nc``,
anywhere below the ROMEX folder (e.g. ``romex_coseqrmx_100m_v1a-003.nc``,
``romex_spirermxeum_274-304_100m_v1a.nc``). Every file is searched for each
day. The files have no global attributes, so who processed them is not
recorded in them. Variables (checked on eepc-0184, 2026-11-27):
``t`` (profile; s since 1970-01-01; 0 in a few profiles, which have no
time), ``altitude`` (500 levels, 10 to 59.9 km every 100 m),
``lat_full``/``lon_full`` (profile x level, degrees), ``lat``/``lon`` (one
per profile), ``temp`` (degrees C), ``press`` (mb), ``ref`` (N),
``Bend_ang`` (rad), ``azim`` (degree), ``is_bad`` (profile; "Is Profile
Bad", 0 or 1; meaning not documented in the files), ``filestamp``
(occultation id, e.g. ``C2E1.2022.244.05.02.R24``).

Output: ``Temp`` = temp + 273.15 (K); ``Pres`` = press (mb = hPa);
``Alt`` = altitude (km): geometric height (not stated in the files; on
2022-10-15 the height increments exceed the hypsometric thickness by the
geometric/geopotential ratio, 0.5% at 10-20 km to 1.5% at 40-50 km, and as
geopotential height COSMIC-2 profiles sit within 30-70 m of ERA5's at the
same pressure); positions per level from
lat_full/lon_full, else the profile's lat/lon. Profiles without a time
(``t`` = 0) are not returned; the number is recorded in the time
variable's comment. With ``instrument_qc``, every level of a profile with
``is_bad`` 1 is ``instrument_flag`` (removed unless ``keep_outliers``).
Other numeric file variables (``ref``, ``Bend_ang``, ``azim``, ``is_bad``,
``lat``, ``lon``) are available through ``additional_vars``. With
``file_source``, ``original_files`` names each profile's file (mission).

ROMEX is kept apart from "GNSS" (the merged daily files), which holds many
of the same occultations for these months.

Reading: for each file, the day's profiles are found from ``t``; they lie
close together in the files (2026-11-27: the block read was at most 1.7 times
the profiles wanted), so one contiguous block of rows is read per file.

Origin: new code (no MATLAB equivalent).

Change log:
2026/11/27 [Claude] initial version
2026/11/27 [Claude] altitude found to be geometric (hypsometric and ERA5 checks)
"""

import math
import os
import re

import numpy as np

from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs._common import find_in, note, note_sources, note_version, record
from atmosgw.curate.obs._limb_tools import conform_vars
from atmosgw.curate.profiles import concat_profiles
from atmosgw.curate.readers import read_netcdf

KELVIN = 273.15
_VERSION = re.compile(r"_(v\w+?)(?:-\d+)?\.nc$")
FILE_NAMES = {"Temp": "temp", "Pres": "press", "Lat": "lat_full", "Lon": "lon_full",
              "Alt": "altitude", "Time": "t"}


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    note_sources(inst_info, dict(FILE_NAMES))
    note(inst_info, "Temp", "file temp (degrees C) + 273.15")
    note(inst_info, "Lat", "lat_full, else the profile's lat")
    note(inst_info, "Lon", "lon_full, else the profile's lon")
    first = math.floor(min(settings["time_range"]))
    last = math.floor(max(settings["time_range"]))
    epoch = datenum(1970, 1, 1)
    span = ((first - epoch) * 86400.0, (last + 1 - epoch) * 86400.0)
    files = sorted(p for p in find_in(inst_info, "romex_*.nc") if os.path.isfile(p))
    if not files:
        record(inst_info)["files_not_found"].append("any romex_*.nc file")
    parts, file_list, no_time = [], [], 0
    for k, path in enumerate(files, start=1):
        if settings["verbose"]:
            print(f"----> Loading ROMEX data: file {k} of {len(files)}")
        seconds = read_netcdf(path, ["t"])["t"]
        no_time += int(np.sum((seconds == 0) & np.isfinite(seconds)))
        with np.errstate(invalid="ignore"):
            rows = np.flatnonzero((seconds > 0) & (seconds >= span[0])
                                  & (seconds < span[1]))
        if rows.size == 0:
            continue
        file_list.append(os.path.basename(path))
        found = _VERSION.search(os.path.basename(path))
        note_version(inst_info, path, f"ROMEX {found.group(1)}" if found else "")
        parts.append(read_rows(path, rows, seconds, var_names, len(file_list),
                               settings))
    if no_time:
        note(inst_info, "Time", f"file t (s since 1970-01-01); profiles with t = 0 "
                                f"are not returned ({no_time} in the files read)")
    else:
        note(inst_info, "Time", "file t (s since 1970-01-01)")
    return concat_profiles(parts, var_names, consume=True), file_list


def read_rows(path, rows, seconds, var_names, file_count, settings):
    """The given profile rows (0-based, sorted) of one file, as [profiles x
    levels] arrays named as in var_names."""
    wanted = {"temp", "press", "lat_full", "lon_full", "lat", "lon"}
    if "InstrumentBad" in var_names:
        wanted.add("is_bad")
    extras = [n for n in settings["additional_vars"]]
    lo, hi = int(rows[0]), int(rows[-1]) + 1
    block = read_netcdf(path, sorted(wanted | set(extras)), rows=slice(lo, hi))
    take = rows - lo
    values = {name: value[take] for name, value in block.items()}
    alt = read_netcdf(path, ["altitude"])["altitude"]
    n_prof, n_lev = rows.size, alt.size
    lat, lon = values["lat_full"], values["lon_full"]
    lat = np.where(np.isnan(lat), values["lat"][:, None], lat)
    lon = np.where(np.isnan(lon), values["lon"][:, None], lon)
    store = {}
    for name in var_names:
        if name == "Temp":
            value = values["temp"] + KELVIN
        elif name == "Pres":
            value = values["press"]
        elif name == "Lat":
            value = lat
        elif name == "Lon":
            value = lon
        elif name == "Alt":
            value = np.tile(alt, (n_prof, 1))
        elif name == "Time":
            value = np.tile((datenum(1970, 1, 1) + seconds[rows] / 86400.0)[:, None],
                            (1, n_lev))
        elif name == "SourceProf":
            value = np.tile((rows + 1.0)[:, None], (1, n_lev))  # 1-based file row
        elif name == "SourceFile":
            value = np.full((n_prof, n_lev), float(file_count))
        elif name == "InstrumentBad":
            value = np.tile((values["is_bad"] == 1).astype(float)[:, None], (1, n_lev))
        elif name in values:
            value = values[name]
        else:
            raise KeyError(f"Variable {name} not found in {path}")
        store[name] = value
    return conform_vars(store, settings)
