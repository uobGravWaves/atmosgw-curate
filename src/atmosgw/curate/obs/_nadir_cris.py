"""
Nadir-sounder reader for CrIS 2-D brightness temperatures (L. Hoffmann's
daily perturbation files), for SNPP and NOAA-20.

Files, below the CrIS folder (``datasets.toml``, ``CrIS``):
``pert_cris_<sat>_<yyyy>_<ddd>.nc`` anywhere below it (on eepc-0184
``2016/`` for SNPP, ``2019/SNPP/`` and ``2019/J1/``), ``<sat>`` "snpp" or
"j1" (NOAA-20), one per day: dimensions ``NTRACK`` (scans: usually 10,800
a day, every 8 s, the first with a mean time of 00:00:12.88 in all 90 files
on eepc-0184; SNPP 2016-01-17 has 8,820, with a 4.6-hour gap),
``NXTRACK`` (30 fields of regard, FOR, across
the scan) and ``NFOV`` (the 9 fields of view, FOV, of each FOR); per
footprint ``time`` (seconds since 2000-01-01T00:00Z; the same for the 9 FOVs
of a FOR), ``lat``, ``lon`` and, for each channel group, ``bt_<group>``
(K) and, for 4.3 and 15 micron, the producer's perturbation
``bt_<group>_pt`` (K) and variance ``bt_<group>_var`` (K^2). Missing values
are NaN, or netCDF's default fill 9.97e36 (no ``_FillValue`` attribute;
SNPP 2016-01-17 and NOAA-20 2019-01-13, where every 15 micron and
perturbation value is fill): any value above 1e30 in size is treated as
-9999 (missing, counted with the fill values). FORs 1 and 30, the
outermost, are NaN in every channel in the files checked (2016-01-01;
2019-01-01 and 15 of both satellites): they are kept on the native grid as
``missing``, and the regridded swath spans FORs 2-29.

Channel groups (``channels`` option; default all, in this order):

========== ===================== ============================= ===========
name       file variable         what (from the files)         nominal
                                                                cm-1
========== ===================== ============================= ===========
4mu        bt_4mu                4.3 micron                    2325.6
8mu        bt_8mu                8.1 micron (no perturbation)  1234.6
10mu       bt_10mu               10.4 micron (no perturbation)  961.5
15mu_low   bt_15mu_low           15 micron, low altitudes       666.7
15mu_high  bt_15mu_high          15 micron, high altitudes      666.7
========== ===================== ============================= ===========

The wavenumbers are nominal (10^4 / wavelength in micron, from the files'
descriptions), not the channels actually averaged: UNVERIFIED.

Layout (the native grid, ``original_grid=True``): the index layout, as the
group's 2022 MATLAB granules: across, 3 x FOR + column; along, 3 x scan +
row, with FOV k (1-based) at row (k - 1) // 3 and column (k - 1) % 3. At
nadir this is a regular 3 x 3 grid about 16 km apart (FOV 1-3 behind the
FOR centre, 7-9 ahead; 1, 4, 7 on the side of FOR 1). Towards the swath
edges the FOV pattern of a FOR rotates (about 45-55 degrees in the files)
and spreads (FOR centres 135 km apart), so neighbours in the array are no
longer along and across the track there (checked 2026-11-27 on the files of
2016-01-01 and 2019-01-01; that the rotation follows the scan angle is the
instrument's design: UNVERIFIED). Each footprint keeps its own position, so
the regridded swath is made from the true positions instead
(``swath.regrid_footprints``).

Granules: pseudo-granules of 6 minutes (45 scans, 135 rows of the index
layout, as AIRS), 240 a day: granule g holds the scans whose mean time is in
the g-th 6-minute slot from the day's nominal first scan (``granule_slots``),
so a granule number means the same time of day in every file, and a file
with a gap gives partial or missing granules (recorded in
``files_not_found`` where the request covers them). Chosen by the mean of
their times, or given explicitly (``granules``). ``satellite`` chooses the
files ("SNPP", the default, or "NOAA-20"; "J1" also accepted).

Variables (reader contract, ``docs/design/nadir_loader.md`` section 5):
``Lat``, ``Lon``, ``Time`` (datenums), ``Temp`` (brightness temperature)
[cross x along x channel], ``Channel`` (wavenumbers), ``ChannelName``, per
row ``Granule``, ``Day``, ``SourceFile``, ``SourceRow`` (1-based row of the
index layout within its granule), and the producer's perturbation and
variance ``BTPert``, ``BTVar`` (NaN for 8.1 and 10.4 micron), declared in
``DerivedVars`` as for AIRS 2-D.

Origin: new code (no MATLAB equivalent).

Change log:
2026/11/27 [Claude] initial version
2026/11/27 [Claude] review: netCDF default fill, granules by time slot, no previous day
"""

import math
import os
import warnings

import netCDF4
import numpy as np

from atmosgw.curate.obs._common import (
    FILL_VALUES,
    day_label,
    locate,
    nadir_fill_to_nan,
    note,
    note_fill_counts,
    note_sources,
    note_version,
    record,
    year_and_day,
)
from atmosgw.curate.obs._nadir_airs import J2000
from atmosgw.curate.readers import read_netcdf

SCANS_PER_GRANULE = 45
GRANULES_PER_DAY = 240
SCAN_S = 8.0  # seconds between scans
FIRST_SCAN_S = 12.88  # mean time of a day's first scan, seconds after midnight
BIG_FILL = 1e30  # |value| above: netCDF's default fill (9.97e36), missing
GRANULE_S = SCANS_PER_GRANULE * SCAN_S  # 360 s
N_FOR, N_FOV = 30, 9
SATELLITES = {"SNPP": "snpp", "NOAA-20": "j1", "J1": "j1"}

# name: (file variable, nominal wavenumber cm-1, description, has perturbation)
CHANNELS = {
    "4mu": ("bt_4mu", 1e4 / 4.3, "4.3 micron", True),
    "8mu": ("bt_8mu", 1e4 / 8.1, "8.1 micron", False),
    "10mu": ("bt_10mu", 1e4 / 10.4, "10.4 micron", False),
    "15mu_low": ("bt_15mu_low", 1e4 / 15, "15 micron, low altitudes", True),
    "15mu_high": ("bt_15mu_high", 1e4 / 15, "15 micron, high altitudes", True),
}

DERIVED = {
    "BTPert": {"output": "brightness_temperature_perturbation", "suffix": "_pt",
               "units": "K",
               "long_name": "brightness temperature perturbation",
               "comment": "from the files: computed by their producer (L. "
                          "Hoffmann's processing); method not checked"},
    "BTVar": {"output": "brightness_temperature_variance", "suffix": "_var",
              "units": "K2",
              "long_name": "brightness temperature variance",
              "comment": "from the files: computed by their producer (L. "
                         "Hoffmann's processing); method not checked"},
}


def satellite_code(satellite):
    """The file-name code ("snpp", "j1") of a satellite name (any case)."""
    name = "SNPP" if satellite is None else str(satellite).upper()
    if name not in SATELLITES:
        raise ValueError(f"satellite must be one of {', '.join(SATELLITES)}, not "
                         f"{satellite!r}.")
    return SATELLITES[name]


def find_file(inst_info, day, code):
    """Path of a day's file anywhere below the folder, or None."""
    year, doy = year_and_day(day)
    name = f"pert_cris_{code}_{year:04d}_{doy:03d}.nc"
    root = inst_info["path"]
    return locate(inst_info, os.path.join(root, name), name, day_label(day))


def index_layout(values):
    """[scan x FOR x FOV (x ...)] -> [cross x along (x ...)] in the index
    layout: cross = 3 FOR + column, along = 3 scan + row, FOV k (0-based) at
    row k // 3, column k % 3."""
    n_scan = values.shape[0]
    out = np.empty((3 * N_FOR, 3 * n_scan) + values.shape[3:], dtype=values.dtype)
    for k in range(N_FOV):
        row, col = divmod(k, 3)
        out[col::3, row::3] = np.swapaxes(values[:, :, k], 0, 1)
    return out


def granule_slots(scan_time, day):
    """Granule number (1-240) of each scan from its mean time (datenums):
    the 6-minute slot it falls in, counted from the day's nominal first scan
    (0 for a scan without a time)."""
    with np.errstate(invalid="ignore"):
        slot = np.floor(((scan_time - day) * 86400 - FIRST_SCAN_S + SCAN_S / 2)
                        / GRANULE_S) + 1
    slot = np.where(np.isfinite(slot), np.clip(slot, 1, GRANULES_PER_DAY), 0)
    return slot.astype(int)


def load(settings, inst_info, var_names):
    """Read the granules; returns (data, file_list). See the module docstring."""
    channels = list(settings.get("channels") or CHANNELS)
    unknown = [c for c in channels if c not in CHANNELS]
    if unknown:
        raise ValueError(f"Unknown CrIS channel(s) {unknown}; choose from "
                         f"{', '.join(CHANNELS)}.")
    code = satellite_code(settings.get("satellite"))
    note_sources(inst_info, {"Lat": "lat", "Lon": "lon", "Time": "time",
                             "Temp": ", ".join(CHANNELS[c][0] for c in channels)})
    note(inst_info, "Time", "time: seconds since 2000-01-01T00:00Z")
    rec = record(inst_info)
    explicit = settings.get("granules") is not None
    if explicit:
        wanted = {}
        for day, g in settings["granules"]:
            wanted.setdefault(int(day), []).append(int(g))
        days = sorted(wanted)
        t0, t1 = -np.inf, np.inf
    else:
        t0, t1 = settings["time_range"]
        # (granule 240 of a day ends at about 00:00:05 the next day, but is
        # chosen by its mean time, about 23:57: the day before never adds one)
        days = list(range(int(math.floor(t0)), int(math.floor(t1)) + 1))
        r0, r1 = settings.get("request_range", settings["time_range"])
    nominal = (FIRST_SCAN_S - SCAN_S / 2
               + (np.arange(1, GRANULES_PER_DAY + 1) - 0.5) * GRANULE_S) / 86400

    granules, file_list, fill_counts = [], [], {}
    for day in days:
        path = find_file(inst_info, day, code)
        if path is None:
            if explicit:
                rec["files_not_found"].extend(f"{day_label(day)} granule {g}"
                                              for g in wanted[day])
            elif np.any((day + nominal >= r0) & (day + nominal <= r1)):
                rec["files_not_found"].append(day_label(day))
            continue
        rec["files_opened"] = rec.get("files_opened", 0) + 1
        seconds = read_netcdf(path, ["time"])["time"]
        with np.errstate(invalid="ignore"):
            valid = np.where((seconds > 0) & (np.abs(seconds) < BIG_FILL), seconds,
                             np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # (scans without times)
            scan_time = J2000 + np.nanmean(valid, axis=(1, 2)) / 86400  # (datenums)
        slot = granule_slots(scan_time, day)
        chosen = []
        for g in range(1, GRANULES_PER_DAY + 1):
            scans = np.flatnonzero(slot == g)
            if explicit:
                if g in wanted[day] and scans.size:
                    chosen.append((g, scans, float(np.mean(scan_time[scans]))))
                continue
            if scans.size:
                mean = float(np.mean(scan_time[scans]))
                if t0 <= mean <= t1:
                    chosen.append((g, scans, mean))
            elif r0 <= day + nominal[g - 1] <= r1:
                rec["files_not_found"].append(f"{day_label(day)} granule {g} "
                                              f"(no scans in {os.path.basename(path)})")
        if explicit:
            for g in sorted(set(wanted[day]) - {c[0] for c in chosen}):
                rec["files_not_found"].append(f"{day_label(day)} granule {g} "
                                              f"(no scans in {os.path.basename(path)})")
        if not chosen:
            continue
        file_list.append(os.path.basename(path))
        note_version(inst_info, path, "")
        first, last = int(chosen[0][1][0]), int(chosen[-1][1][-1]) + 1
        names = ["time", "lat", "lon"]
        for c in channels:
            base = CHANNELS[c][0]
            names.append(base)
            if CHANNELS[c][3]:
                names += [base + d["suffix"] for d in DERIVED.values()]
        values = _read_present(path, names, slice(first, last))
        for g, scans, mean in chosen:
            take = scans - first
            part, counts = _granule(values, take, day, g, channels,
                                    len(file_list))
            for key, n in counts.items():
                fill_counts[key] = fill_counts.get(key, 0) + n
            granules.append((mean, part))

    note_fill_counts(inst_info, fill_counts)
    granules.sort(key=lambda item: item[0])
    names = ["Lat", "Lon", "Time", "Temp", "Granule", "Day", "SourceFile",
             "SourceRow"] + list(DERIVED)
    if granules:
        data = {name: np.concatenate([p[name] for _, p in granules],
                                     axis=1 if granules[0][1][name].ndim >= 2 else 0)
                for name in names}
    else:
        n_x = 3 * N_FOR
        data = {name: np.empty((n_x, 0, len(channels))) if name in ("Temp", *DERIVED)
                else (np.empty((n_x, 0)) if name in ("Lat", "Lon", "Time")
                      else np.empty(0)) for name in names}
    data["Channel"] = np.array([CHANNELS[c][1] for c in channels])
    data["ChannelName"] = np.array(channels)
    data["DerivedVars"] = DERIVED
    plain = [c for c in channels if not CHANNELS[c][3]]
    if plain:
        note(inst_info, "BTPert", "no perturbation or variance in the files for "
                                  f"{', '.join(plain)} (NaN)")
    return data, file_list


def _granule(values, take, day, g, channels, file_number):
    """One pseudo-granule in the index layout: (part, fill counts)."""
    sub = {name: np.array(v[take], dtype=float) for name, v in values.items()}
    for v in sub.values():  # (netCDF's default fill: counted as the others)
        with np.errstate(invalid="ignore"):
            v[np.abs(v) > BIG_FILL] = FILL_VALUES[-1]
    lat, lon = index_layout(sub["lat"]), index_layout(sub["lon"])
    seconds = index_layout(sub["time"])
    # fill values and bad times, as the other nadir readers (which take
    # [along x cross]: transposed views, changed in place)
    span = ((day - J2000) * 86400, (day + 2 - J2000) * 86400)
    counts = nadir_fill_to_nan(lat.T, lon.T, seconds.T,
                               np.zeros(lat.T.shape + (1,)), span)
    planes = {"": [], **{d["suffix"]: [] for d in DERIVED.values()}}
    for c in channels:
        base = CHANNELS[c][0]
        for suffix in planes:
            name = base + suffix
            if name in sub:
                v = index_layout(sub[name])
                fill = np.isin(v, FILL_VALUES)
                counts["temp"] += int(fill.sum())
                v[fill] = np.nan
            else:
                v = np.full(lat.shape, np.nan)
            planes[suffix].append(v)
    n_along = lat.shape[1]
    part = {"Lat": lat, "Lon": lon, "Time": J2000 + seconds / 86400,
            "Temp": np.stack(planes[""], axis=-1),
            "Granule": np.full(n_along, float(g)), "Day": np.full(n_along, float(day)),
            "SourceFile": np.full(n_along, float(file_number)),
            "SourceRow": np.arange(1, n_along + 1, dtype=float)}
    for key, d in DERIVED.items():
        part[key] = np.stack(planes[d["suffix"]], axis=-1)
    return part, counts


def _read_present(path, names, rows):
    """The variables of names that the file has (others left out), rows of
    the first dimension only."""
    with netCDF4.Dataset(path) as f:
        present = [n for n in names if n in f.variables]
    return read_netcdf(path, present, rows=rows)
