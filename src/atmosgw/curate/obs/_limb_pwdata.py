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
``*<identifier>*<yyyy>d<ddd>*``.

The MATLAB module concatenates days in chunks of 15 for speed; the result is
the same as concatenating all days at once, which is done here. One detail is
kept: when exactly one file is read, ``Temp_PW`` keeps the number type it has
in the file (the MATLAB code copies it directly rather than through
``cat_struct``, which converts to double). If a file held single-precision
``Temp_PW``, MATLAB would then interpolate it in single precision, which is
not reproduced here; the HIRDLS and SABER PW files checked (2007d001,
2017d050) hold double precision, so this does not arise for them.
"""

import os

import numpy as np

from atmosgw.curate.atmosphere import height_to_pressure
from atmosgw.curate.obs._limb_tools import conform_vars, year_and_day
from atmosgw.curate.profiles import concat_profiles
from atmosgw.curate.readers import find_files, read_mat


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    first, last = (int(np.floor(t)) for t in settings["time_range"])
    days, file_list = [], []
    for k, day in enumerate(range(first, last + 1), start=1):
        if settings["verbose"]:
            print(
                f"----> Loading Hindley23-style PW data: day {k} of {last - first + 1}"
            )
        year, doy = year_and_day(day)
        if inst_info["inst"].lower() != "misc":
            path = os.path.join(
                inst_info["path"],
                f"{year:04d}",
                f"{inst_info['inst']}_PWs_{year:04d}d{doy:03d}.mat",
            )
        else:
            found = find_files(
                settings["misc_path"],
                f"*{settings['misc_identifier']}*{year:04d}d{doy:03d}*",
            )
            if not found:
                continue
            path = found[0]
        if not os.path.isfile(path):
            continue

        inst_data = read_mat(path)
        file_list.append(os.path.basename(path))
        file_count = len(file_list)

        store = {}
        for name in var_names:
            if name == "Pres" and "Pres" not in inst_data:
                inst_data["Pres"] = height_to_pressure(inst_data["Alt"])
            if name == "SourceProf":
                n_prof, n_lev = np.shape(inst_data["Lat"])[:2]
                profile = np.arange(1, n_prof + 1, dtype=float)  # 1-based
                inst_data["SourceProf"] = np.tile(profile.reshape(-1, 1), (1, n_lev))
            if name == "SourceFile":
                inst_data["SourceFile"] = np.full(np.shape(inst_data["Lat"]),
                                                  float(file_count))
            if name not in inst_data:
                raise KeyError(f"Variable {name} not found in {path}")
            store[name] = inst_data[name]
        days.append(conform_vars(store, settings))

    data = concat_profiles(days, var_names)
    if len(days) == 1:
        data["Temp_PW"] = days[0]["Temp_PW"]  # keeps the file's number type
    return data, file_list
