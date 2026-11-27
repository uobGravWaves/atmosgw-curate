"""
Limb-sounder reader for LIMS (Nimbus-7 Limb Infrared Monitor of the
Stratosphere), Version 6 revised Level 2 ASCII files.

Files: one per day, ``LIMS_2_<yyyymmdd>_V6.ASCII``, anywhere below the LIMS
folder. Format (LIMS V6 User's Guide, section 5, "Data Format And Example";
checked against the files on 2026-11-27): a description block (``####  Year =
1978`` ...), then one block per scan: three marker lines (``*###  SCAN n``),
a header line, six channel lines, and ``nl_std`` (109) lines of 18 values on
the standard pressure grid (1000 to 0.001 mb), top level first.

- Header line: fixed width (195 characters). Fields run together when wide
  (e.g. ``gt2`` and ``icloud`` as ``0.00100``) and the time has blanks
  (``0: 2: 8``), so it is read by column (``HEADER_COLUMNS``), not by
  splitting on blanks. The columns were found from the files: every header
  of the 201 files held on eepc-0184 parses to values in range.
- Level lines: ``dang_std zt_std pt_std tt_std rad_std1..6 tg_std1..2
  qmix_std1..5 geo_heights`` (``LEVEL_NAMES``); 1e24 is missing (NaN).

Output: ``Temp`` = tt_std (K), ``Alt`` = zt_std (km, "height"), ``Pres`` =
pt_std (mb = hPa), position = the 30 km tangent point (``alat``, ``alon``,
degrees east, 0-360) for every level, as the files give one per scan.
``Time`` = the day ``iday`` (day of the year, of the file's year, or the next
year after 31 December) + ``time`` (GMT). Geopotential height
(``geo_heights``, km, referenced to the NMC 50 mb heights) and every other
level or header value are available through ``additional_vars``, by the
names in ``LEVEL_NAMES`` and ``HEADER_COLUMNS``.

Days: the files hold complete orbits, so a file usually starts at the first
orbit of its day and continues into the next; but two start before their
day (1979-03-27 by 1.24 h, 1979-05-15 by 0.12 h; the files do not
overlap). Each request also reads the previous and the next day's files,
and the core keeps the profiles whose mean time is in the requested range.

No screening: the User's Guide says the temperature retrievals are not
screened (section 2) and not cloud-screened (section 3).

Origin: new code (no MATLAB equivalent).

Change log:
2026/11/27 [Claude] initial version
2026/11/27 [Claude] next day's file also read (two files start before their day)
"""

import math
import re

import numpy as np

from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs._common import (
    choose_file,
    day_label,
    find_in,
    note,
    note_sources,
    note_version,
    record,
)
from atmosgw.curate.obs._limb_tools import conform_vars, loop_days

FILL = 1e24  # missing value in the files (User's Guide, section 5 example)
N_CHANNELS = 6  # channel lines after each header (CO2N, CO2W, O3, HNO3, H2O, NO2)

# values on each level line, in order (User's Guide, section 5)
LEVEL_NAMES = ("dang_std", "zt_std", "pt_std", "tt_std",
               "rad_std1", "rad_std2", "rad_std3", "rad_std4", "rad_std5", "rad_std6",
               "tg_std1", "tg_std2", "qmix_std1", "qmix_std2", "qmix_std3",
               "qmix_std4", "qmix_std5", "geo_heights")

# header line: name -> (first, last + 1) character (0-based, half-open), for the
# fields read; the rest of the line is split on blanks (REST_NAMES)
HEADER_COLUMNS = {"nl_std": (0, 5), "ngs1": (5, 9), "nch": (9, 13),
                  "alat": (13, 22), "alon": (22, 31), "iorbit": (31, 36),
                  "irec": (36, 40), "iday": (40, 44), "hour": (44, 47),
                  "minute": (48, 50), "second": (51, 53), "idn": (53, 59),
                  "iud": (59, 62), "alt": (62, 72), "gt1": (72, 78),
                  "gt2": (78, 84), "icloud": (84, 87), "iad": (87, 90)}
REST_NAMES = ("sunasc", "sundec", "grncha", "avrl", "tkm", "plat", "plon", "szad",
              "shift")

# output names -> LIMS names
FILE_NAMES = {"Lat": "alat", "Lon": "alon", "Alt": "zt_std", "Pres": "pt_std",
              "Temp": "tt_std", "Time": "iday, time"}

_YEAR = re.compile(r"^####\s+Year\s*=\s*(\d+)", re.MULTILINE)
_VERSION = re.compile(r"^####\s+Version\s*=\s*(\S+)", re.MULTILINE)


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    note_sources(inst_info, dict(FILE_NAMES))
    note(inst_info, "Time", "file iday (day of the year) + time (GMT)")
    note(inst_info, "Lat", "30 km tangent point, used for the whole profile")
    note(inst_info, "Lon", "30 km tangent point, used for the whole profile")
    note(inst_info, "Alt", "zt_std, height (km); geopotential height is the "
                           "additional variable geo_heights")
    first = math.floor(min(settings["time_range"]))
    last = math.floor(max(settings["time_range"]))

    def find_file(day):
        stamp = day_label(day).replace("-", "")  # yyyymmdd
        found = find_in(inst_info, f"LIMS_2_{stamp}_*")
        path = choose_file(found, inst_info, day_label(day))
        if path is None and first <= day <= last:  # (the days either side help)
            record(inst_info)["files_not_found"].append(day_label(day))
        return path

    def read_file(path, file_count, day):
        scans = read_lims_file(path)
        note_version(inst_info, path, scans["version"])
        n_prof, n_lev = scans["n_scans"], scans["n_levels"]
        store = {}
        for name in var_names:
            if name == "Time":
                value = scans["time"]
            elif name == "SourceProf":
                value = np.arange(1.0, n_prof + 1)  # 1-based scan number
            elif name == "SourceFile":
                value = np.full(n_prof, float(file_count))
            elif name in FILE_NAMES and FILE_NAMES[name] in scans:
                value = scans[FILE_NAMES[name]]
            elif name in scans:
                value = scans[name]
            else:
                raise KeyError(f"Variable {name} not found in {path}")
            # one value per scan: repeated over levels (extras: by conform_vars)
            if np.ndim(value) == 1 and name not in settings["additional_vars"]:
                value = np.tile(value[:, None], (1, n_lev))
            store[name] = value
        return conform_vars(store, settings)

    # from the day before the start to the day after the end (files run into
    # the next day, and some start before their own); a missing file is
    # recorded by find_file, for requested days only
    shifted = dict(settings, time_range=(min(settings["time_range"]) - 1,
                                         max(settings["time_range"]) + 1))
    return loop_days(shifted, var_names, find_file, read_file, "LIMS")


def read_lims_file(path):
    """
    Every scan of one LIMS V6 ASCII file. Returns a dict: ``n_scans``,
    ``n_levels``; each level value (``LEVEL_NAMES``) as [scans x levels]
    float64, 1e24 as NaN; each header value (``HEADER_COLUMNS``,
    ``REST_NAMES``) as a vector over scans, float64; ``time`` (datenums);
    ``version`` (the file's ``Version`` line).
    """
    with open(path, encoding="ascii", errors="replace") as handle:
        text = handle.read()
    year = int(_YEAR.search(text).group(1))
    found = _VERSION.search(text)
    version = f"Version {found.group(1)}" if found else ""
    lines = text.splitlines()
    starts = [k for k, line in enumerate(lines) if line.startswith("*###  SCAN")]
    n_scans = len(starts)
    header = {name: np.full(n_scans, np.nan) for name in
              list(HEADER_COLUMNS) + list(REST_NAMES)}
    level_text = []
    n_levels = None
    for j, k in enumerate(starts):
        line = lines[k + 2]  # (after the marker line below "SCAN n")
        for name, (a, b) in HEADER_COLUMNS.items():
            header[name][j] = float(line[a:b])
        rest = line[HEADER_COLUMNS["iad"][1]:].split()
        for name, value in zip(REST_NAMES, rest):
            header[name][j] = float(value)
        n = int(header["nl_std"][j])
        if n_levels is None:
            n_levels = n
        elif n != n_levels:
            raise ValueError(f"{path}: scan {j + 1} has {n} levels, not {n_levels}")
        first = k + 3 + N_CHANNELS
        level_text.extend(lines[first:first + n])
    n_levels = n_levels or 0
    values = np.array(" ".join(level_text).split(), dtype=float)
    if values.size != n_scans * n_levels * len(LEVEL_NAMES):
        raise ValueError(f"{path}: {values.size} level values, expected "
                         f"{n_scans} x {n_levels} x {len(LEVEL_NAMES)}")
    values = values.reshape(n_scans, n_levels, len(LEVEL_NAMES))
    values[values >= FILL] = np.nan
    out = {name: values[:, :, i] for i, name in enumerate(LEVEL_NAMES)}
    out.update(header)
    out["n_scans"], out["n_levels"], out["version"] = n_scans, n_levels, version
    # day of the year of the file's year; days after 31 December continue
    # into the next year (a file runs into the next day)
    day = np.asarray(header["iday"])
    file_day = day[0] if n_scans else 0
    year_length = datenum(year + 1, 1, 1) - datenum(year, 1, 1)
    day = np.where(day < file_day - 300, day + year_length, day)
    out["time"] = (datenum(year, 1, 1) + day - 1 + header["hour"] / 24
                   + header["minute"] / 1440 + header["second"] / 86400)
    return out
