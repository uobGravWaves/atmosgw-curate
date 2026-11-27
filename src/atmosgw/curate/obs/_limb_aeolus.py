"""
Limb-loader reader for Aeolus (ALADIN) Level 2B horizontal line-of-sight
(HLOS) winds, from monthly "flat" netCDF files made from the ESA L2B product
(ALD_U_N_2B) through VirES by the group's ``get_aeolus.py`` (held beside the
data).

Files: one per month, ``Aeolus_L2B_monthly_<yyyy>_<mm>.nc``, anywhere below
the Aeolus folder. Variables, for each channel (``ray`` Rayleigh, ``mie``
Mie), one value per observation (range bin): ``<ch>_time`` (s since
2000-01-01), ``<ch>_lat``, ``<ch>_lon``, ``<ch>_z_bot``, ``<ch>_z_top`` (m),
``<ch>_hlos``, ``<ch>_hlos_err`` (cm/s; 32767 and 65535 are missing),
``<ch>_obs_type`` (1 cloudy, 2 clear: UNVERIFIED, inferred from the files:
Rayleigh type 2 observations form full columns, the Mie file holds only type
1), ``<ch>_valid`` (0 or 1), ``<ch>_az`` (degrees). The global attribute
``Sources`` lists each L2B file with its baseline and processor version.

Channel (``aeolus_channel``): "rayleigh" (default) uses Rayleigh-clear
observations, "mie" Mie-cloudy (the usual choice for L2B winds; the other
types are not read).

Profiles: the flat files carry no profile identifier, so profiles are
rebuilt (a heuristic, new code): the observations of one channel and type,
in time order (stable), start a new profile where (a) the range bin is more
than 3 km above the previous one (the column starting again at its top),
(b) a range bin (the same bottom and top) is already in the profile, or (c)
12 s (about one basic repeat cycle, 87 km) have passed since the profile's
first observation. On 2019-01 and 2021-06 Rayleigh-clear data, rule (a)
alone gave profiles spanning at most 11.6 s. Within a profile the levels
are ordered by height (bottom first). Each level keeps its own time and
position. The rebuild is approximate: on 2019-01-05 (Rayleigh-clear) it made
about 2.2 profiles per 12 s rather than one per cycle, and 24% of the
observations (median height 1.9 km) are in profiles of 10 levels or fewer,
because near-surface bins have times of their own, before or after their
column's (independent review, 2026-11-27). A profile identifier from the L2B
product would make it exact.

Duplicates: consecutive L2B files overlap by about 30 s, so the flat files
hold copies of some observations; an observation identical in time, bin,
type, HLOS and error to another is read once (2019-01-05, Rayleigh-clear:
758 of 164,750 share time and bin with another, 180 of them also the wind; the
other pairs differ in wind and are both kept). The
observations read reach 12 s beyond the requested days, so a profile
crossing midnight is whole (the core selects profiles by mean time); one
crossing the start of a file (about 00:00:35 on the first day of a file) is
split in two.

Output: ``HLOS`` (m s-1) and ``HLOSErr`` (m s-1, the L2B error estimate as
distributed) from cm/s; ``LOSAz`` = ``<ch>_az``; ``Alt`` = the middle of the
range bin (km); pressure estimated from height (no pressure in the files).
With ``instrument_qc``, observations with ``valid`` 0 are
``instrument_flag``. Extras (``additional_vars``): ``z_bot``, ``z_top``
(km), ``obs_type``, ``valid``.

Sign convention of HLOS: PLACEHOLDER (as distributed in the L2B product; not
yet checked). The group's MATLAB ``aeolus_uv2hlos.m`` computes
HLOS = -u sin(az) - v cos(az).

Origin: new code (no MATLAB equivalent).

Change log:
2026/11/27 [Claude] initial version
2026/11/27 [Claude] exact duplicate observations read once; 12 s margin at day edges
"""

import math
import re

import numpy as np

from atmosgw.curate.atmosphere import height_to_pressure
from atmosgw.curate.matlab_time import datenum, datenum_to_datetime
from atmosgw.curate.obs._common import (
    choose_file,
    find_in,
    note,
    note_sources,
    note_version,
    record,
)
from atmosgw.curate.obs._limb_tools import conform_vars, note_estimated_pres
from atmosgw.curate.profiles import concat_profiles
from atmosgw.curate.readers import read_attrs, read_netcdf

CHANNELS = {"rayleigh": ("ray", 2), "mie": ("mie", 1)}  # (prefix, obs_type used)
JUMP_KM = 3.0  # a bin this far above the previous one starts a new profile
MAX_SPAN_S = 12.0  # a profile spans at most this long
HLOS_FILL, ERR_FILL = 32767, 65535
SIGN_CONVENTION = ("PLACEHOLDER: as distributed in the Aeolus L2B product; "
                   "not yet checked")
EXTRAS = ("z_bot", "z_top", "obs_type", "valid")

_SOURCE = re.compile(r"'(\w+)',\s*'([\w.]+)'\)")


def load(settings, inst_info, var_names):
    """Read all months; returns (data, file_list). See limb.py for the arguments."""
    channel = str(settings.get("aeolus_channel") or "rayleigh").lower()
    prefix, obs_type = CHANNELS[channel]
    note_sources(inst_info, {"HLOS": f"{prefix}_hlos", "HLOSErr": f"{prefix}_hlos_err",
                             "LOSAz": f"{prefix}_az", "Lat": f"{prefix}_lat",
                             "Lon": f"{prefix}_lon", "Time": f"{prefix}_time",
                             "Alt": f"{prefix}_z_bot, {prefix}_z_top"})
    kind = "clear" if obs_type == 2 else "cloudy"
    note(inst_info, "HLOS", f"{channel.capitalize()}-{kind} observations, cm/s to "
                            f"m s-1; sign convention {SIGN_CONVENTION}")
    note(inst_info, "HLOSErr", "L2B HLOS error estimate as distributed, cm/s to m s-1")
    note(inst_info, "Alt", "middle of the range bin")
    note(inst_info, "Time", "each observation's own time (s since 2000-01-01)")
    note_estimated_pres(inst_info)

    first = math.floor(min(settings["time_range"]))
    last = math.floor(max(settings["time_range"]))
    # the "monthly" files are made of weekly downloads, so they run past the
    # month's end (e.g. 2018_12 to 2019-01-05): each month's days are looked
    # for in its own file and the previous month's
    requested = list(dict.fromkeys(_year_month(d) for d in range(first, last + 1)))
    months = []
    for year, mon in requested:
        previous = (year - 1, 12) if mon == 1 else (year, mon - 1)
        months += [m for m in (previous, (year, mon)) if m not in months]
    parts, file_list = [], []
    for k, (year, mon) in enumerate(months, start=1):
        if settings["verbose"]:
            print(f"----> Loading Aeolus data: file {k} of {len(months)}")
        found = find_in(inst_info, f"Aeolus_L2B_monthly_{year:04d}_{mon:02d}.nc")
        path = choose_file(found, inst_info, f"{year:04d}-{mon:02d}")
        if path is None:
            if (year, mon) in requested:
                record(inst_info)["files_not_found"].append(f"{year:04d}-{mon:02d}")
            continue
        part = read_month(path, prefix, obs_type, (first, last + 1), var_names,
                          len(file_list) + 1, settings)
        if np.shape(part["Lat"])[0] == 0 and (year, mon) not in requested:
            continue  # (the previous month's file holds none of these days)
        file_list.append(path.rsplit("/", 1)[-1])
        note_version(inst_info, path, source_version(path))
        parts.append(part)
    return concat_profiles(parts, var_names, consume=True), file_list


def read_month(path, prefix, obs_type, day_span, var_names, file_count, settings):
    """
    The profiles of one month file whose observations fall in day_span
    (datenums [start, end)), as [profiles x levels] arrays named as in
    var_names (Lat, Lon, Time, HLOS, HLOSErr, LOSAz, Pres, Alt, SourceProf,
    SourceFile, PresEstimated, InstrumentBad and extras).
    """
    names = [f"{prefix}_{n}" for n in ("time", "lat", "lon", "z_bot", "z_top", "hlos",
                                       "hlos_err", "obs_type", "valid", "az")]
    raw = read_netcdf(path, [names[0], names[7]])
    seconds = raw[names[0]]
    epoch = datenum(2000, 1, 1)
    start_s, end_s = ((d - epoch) * 86400.0 for d in day_span)
    start_s, end_s = start_s - MAX_SPAN_S, end_s + MAX_SPAN_S  # (whole profiles)
    with np.errstate(invalid="ignore"):
        rows = np.flatnonzero((raw[names[7]] == obs_type) & (seconds >= start_s)
                              & (seconds < end_s))
    if rows.size:
        rows = rows[np.argsort(seconds[rows], kind="stable")]
        lo, hi = rows.min(), rows.max() + 1  # (one contiguous read)
        block = read_netcdf(path, names, rows=slice(lo, hi))
        obs = {n.split("_", 1)[1]: block[n][rows - lo] for n in names}
    else:
        obs = {n.split("_", 1)[1]: np.empty(0) for n in names}
    if obs["time"].size:  # (exact copies from overlapping L2B files: once)
        key = np.stack([np.nan_to_num(obs[n], nan=-1.0) for n in
                        ("time", "z_bot", "z_top", "obs_type", "hlos", "hlos_err")])
        order = np.lexsort(key[::-1])
        new = np.r_[True, np.any(np.diff(key[:, order], axis=1) != 0, axis=0)]
        keep = np.sort(order[new])  # (the first of each set, in time order)
        obs = {n: v[keep] for n, v in obs.items()}
    obs["hlos"][obs["hlos"] == HLOS_FILL] = np.nan
    obs["hlos_err"][obs["hlos_err"] == ERR_FILL] = np.nan
    group = profile_groups(obs["time"], obs["z_bot"] / 1000, obs["z_top"] / 1000)

    alt = (obs["z_bot"] + obs["z_top"]) / 2000.0  # km
    per_point = {
        "Lat": obs["lat"], "Lon": obs["lon"],
        "Time": epoch + obs["time"] / 86400.0,
        "HLOS": obs["hlos"] / 100.0, "HLOSErr": obs["hlos_err"] / 100.0,
        "LOSAz": obs["az"], "Alt": alt, "Pres": height_to_pressure(alt),
        "PresEstimated": np.ones(alt.shape),
        "InstrumentBad": (obs["valid"] == 0).astype(float),
        "z_bot": obs["z_bot"] / 1000.0, "z_top": obs["z_top"] / 1000.0,
        "obs_type": obs["obs_type"], "valid": obs["valid"],
    }
    n_prof = int(group.max()) + 1 if group.size else 0
    counts = np.bincount(group, minlength=n_prof) if n_prof else np.zeros(0, int)
    n_lev = int(counts.max()) if n_prof else 0
    # level within each profile, bottom first (by height, stable)
    order = np.lexsort((alt, group)) if group.size else np.zeros(0, int)
    starts = np.r_[0, np.cumsum(counts)[:-1]] if n_prof else np.zeros(0, int)
    level = np.arange(order.size) - np.repeat(starts, counts)
    store = {}
    for name in var_names:
        if name == "SourceProf":
            value = np.arange(1.0, n_prof + 1)[group]  # 1-based, in this file
        elif name == "SourceFile":
            value = np.full(group.shape, float(file_count))
        elif name in per_point:
            value = np.asarray(per_point[name], dtype=float)
        else:
            raise KeyError(f"Variable {name} not found in {path}")
        grid = np.full((n_prof, n_lev), np.nan)
        grid[group[order], level] = value[order]
        store[name] = grid
    return conform_vars(store, settings)


def profile_groups(seconds, z_bot, z_top):
    """
    Profile number (0-based) of each observation, given in time order: a new
    profile starts where the bin is more than JUMP_KM above the previous
    one, where the bin (z_bot, z_top) is already in the profile, or where
    MAX_SPAN_S have passed since the profile's first observation.
    """
    n = np.size(seconds)
    group = np.zeros(n, dtype=np.int64)
    current, seen, began = 0, set(), None
    for k in range(n):
        key = (z_bot[k], z_top[k])
        if k > 0 and (z_bot[k] - z_bot[k - 1] > JUMP_KM or key in seen
                      or seconds[k] - began > MAX_SPAN_S):
            current, seen = current + 1, set()
        if not seen:
            began = seconds[k]
        seen.add(key)
        group[k] = current
    return group


def _year_month(day):
    """(year, month) of an integer datenum."""
    stamp = str(datenum_to_datetime(day))
    return int(stamp[:4]), int(stamp[5:7])


def source_version(path):
    """The L2B baselines and processor versions the file's Sources list, e.g.
    "L2B baseline 2B16 processor 03.90"."""
    text = str(read_attrs(path).get("Sources", ""))
    found = sorted(set(_SOURCE.findall(text)))
    return "; ".join(f"L2B baseline {b} processor {p}" for b, p in found)
