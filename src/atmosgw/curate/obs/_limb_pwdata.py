"""
Limb-sounder reader for planetary-wave-filtered data ("Hindley23" format,
made by the group's MATLAB code). Port of ``module_load_pwdata.m``
(uobGravWaves/MultiInstrumentFunctions, commit 78eef78).

File format: one MATLAB .mat file per day with [profiles x levels] arrays
``Lat``, ``Lon``, ``Time``, ``Alt``, ``Temp_Residual`` and ``Temp_PW``
([profiles x levels x modes]); ``Pres`` is estimated from ``Alt`` with
``height_to_pressure`` if absent.

File names: ``<pws path>/<yyyy>/<INSTRUMENT>_PWs_<yyyy>d<ddd>.mat``, or for
Misc the first file below ``misc_path`` matching
``*<identifier>*<yyyy>d<ddd>*``. A file not at its usual path is looked for
by name anywhere below the pws folder (a flat folder, or other nesting).

The MATLAB module concatenates days in chunks of 15 for speed; the result is
the same as concatenating all days at once, which is done here. One detail is
kept: when exactly one file is read, ``Temp_PW`` keeps the number type it has
in the file (the MATLAB code copies it directly rather than through
``cat_struct``, which converts to double). If a file held single-precision
``Temp_PW``, MATLAB would then interpolate it in single precision, which is
not reproduced here; the HIRDLS and SABER PW files checked (2007d001,
2017d050) hold double precision, so this does not arise for them.

GNSS PW files get the same corrections as the raw GNSS files
(``_limb_ace_gnss.gnss_fixes``: fill values, times, positions, repeated
occultations); their pressure is estimated from height, as before.
"""

import os
import re

import numpy as np

from atmosgw.curate.atmosphere import height_to_pressure
from atmosgw.curate.obs._common import (
    choose_file,
    day_label,
    find_in,
    locate,
    note_sources,
    note_version,
    record,
    year_and_day,
)
from atmosgw.curate.obs._limb_ace_gnss import gnss_fixes, note_gnss_counts
from atmosgw.curate.obs._limb_tools import conform_vars, note_estimated_pres
from atmosgw.curate.profiles import concat_profiles
from atmosgw.curate.readers import read_mat


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    first, last = (int(np.floor(t)) for t in settings["time_range"])
    note_sources(inst_info, {name: name for name in var_names
                             if name not in ("SourceProf", "SourceFile")})
    days, file_list = [], []
    for k, day in enumerate(range(first, last + 1), start=1):
        if settings["verbose"]:
            print(
                f"----> Loading Hindley23-style PW data: day {k} of {last - first + 1}"
            )
        year, doy = year_and_day(day)
        if inst_info["inst"].lower() != "misc":
            name = f"{inst_info['inst']}_PWs_{year:04d}d{doy:03d}.mat"
            path = locate(inst_info, os.path.join(inst_info["path"], f"{year:04d}",
                                                  name), name, day_label(day))
        else:
            found = find_in(
                inst_info,
                f"*{settings['misc_identifier']}*{year:04d}d{doy:03d}*",
                root=settings["misc_path"],
            )
            path = choose_file(found, inst_info, day_label(day))
        if path is None or not os.path.isfile(path):
            record(inst_info)["files_not_found"].append(day_label(day))
            continue

        inst_data = read_mat(path)
        file_list.append(os.path.basename(path))
        file_count = len(file_list)
        note_version(inst_info, path, _pw_version(inst_data))
        rows = np.arange(np.shape(inst_data["Lat"])[0])  # file rows (0-based)
        if inst_info["inst"] == "GNSS":  # (as the raw GNSS files)
            rows = gnss_fixes(inst_data, inst_info,
                              settings.get("keep_duplicate_occultations", False))

        store = {}
        estimated = "Pres" not in inst_data
        for name in var_names:
            if name == "Pres" and estimated:
                inst_data["Pres"] = height_to_pressure(inst_data["Alt"])
                note_estimated_pres(inst_info)
            if name == "PresEstimated":  # (for the output's data states)
                inst_data[name] = np.full(np.shape(inst_data["Alt"]), float(estimated))
            if name == "SourceProf":
                n_lev = np.shape(inst_data["Lat"])[1]
                profile = (rows + 1).astype(float)  # 1-based row in the file
                inst_data["SourceProf"] = np.tile(profile.reshape(-1, 1), (1, n_lev))
            if name == "SourceFile":
                inst_data["SourceFile"] = np.full(np.shape(inst_data["Lat"]),
                                                  float(file_count))
            if name not in inst_data:
                raise KeyError(f"Variable {name} not found in {path}")
            store[name] = inst_data[name]
        days.append(conform_vars(store, settings))

    note_gnss_counts(inst_info)
    data = concat_profiles(days, var_names)
    if len(days) == 1:
        data["Temp_PW"] = days[0]["Temp_PW"]  # keeps the file's number type
    return data, file_list


def _pw_version(inst_data):
    """'PW code version <date>' from the file's Note, if it gives one."""
    note_text = " ".join(str(x) for x in np.asarray(inst_data.get("Note", ""),
                                                     dtype=object).ravel())
    found = re.search(r"Generated using (\S+) code version", note_text)
    return f"PW code version {found.group(1)}" if found else ""
