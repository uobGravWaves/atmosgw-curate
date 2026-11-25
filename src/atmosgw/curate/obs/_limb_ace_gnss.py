"""
Limb-sounder reader for ACE-FTS, GNSS radio occultation, and "Misc" files
(any files in the same MATLAB format). Port of ``module_load_ACE_GNSS.m``
(uobGravWaves/MultiInstrumentFunctions, commit 78eef78).

File format: one MATLAB .mat file per day, holding [profiles x levels] arrays
named as in the output (``Lat``, ``Lon``, ``Temp``, and for ACE and Misc
``Alt`` and ``Time``). GNSS files instead hold ``MSL_alt`` (one row of
heights) and a table ``MetaData`` with a ``time`` column (one per profile).
If there is no ``Pres``, it is estimated from ``Alt`` with
``height_to_pressure``.

File names: ACE ``<path>/<yyyy>/ace_v52_<yyyy>d<ddd>.mat``; GNSS
``<path>/<yyyy>/merged_ro_<yyyy>_<ddd>.mat``; Misc: the first file below
``misc_path`` matching ``*<identifier>*<yyyy>d<ddd>*.mat``.
"""

import os

import numpy as np

from atmosgw.curate.atmosphere import height_to_pressure
from atmosgw.curate.obs._common import (
    choose_file,
    day_label,
    note,
    note_sources,
    year_and_day,
)
from atmosgw.curate.obs._limb_tools import ESTIMATED_PRES, conform_vars, loop_days
from atmosgw.curate.readers import find_files, read_mat


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    instrument = settings["instrument"]
    sources = {name: name for name in ("Lat", "Lon", "Temp", "Pres", "Alt", "Time")}
    if instrument == "GNSS":
        sources.update(Alt="MSL_alt", Time="MetaData.time")
    note_sources(inst_info, sources)

    def find_file(day):
        year, doy = year_and_day(day)
        if instrument == "GNSS":
            return os.path.join(
                inst_info["path"], f"{year:04d}", f"merged_ro_{year:04d}_{doy:03d}.mat"
            )
        if instrument == "ACE":
            return os.path.join(
                inst_info["path"], f"{year:04d}", f"ace_v52_{year:04d}d{doy:03d}.mat"
            )
        found = find_files(
            inst_info["path"],
            f"*{inst_info['identifier']}*{year:04d}d{doy:03d}*.mat",
        )
        return choose_file(found, inst_info, day_label(day))

    def read_file(path, file_count, day):
        inst_data = read_mat(path)
        if "Pres" not in inst_data:
            inst_data["Pres"] = height_to_pressure(inst_data["Alt"])
            note(inst_info, "Pres", ESTIMATED_PRES + ", for files without Pres")

        store = {}
        for name in var_names:
            computed = ("Alt", "Time", "SourceProf", "SourceFile")
            if name not in computed and name not in inst_data:
                raise KeyError(f"Variable {name} not found, terminating")
            if name == "Time":
                if instrument == "GNSS":
                    time = np.asarray(inst_data["MetaData"]["time"], dtype=float)
                    n_lev = np.shape(inst_data["Lat"])[1]
                    store["Time"] = np.tile(time.reshape(-1, 1), (1, n_lev))
                else:
                    store["Time"] = inst_data["Time"]
            elif name == "Alt":
                if instrument == "GNSS":
                    n_prof = np.shape(inst_data["Lat"])[0]
                    store["Alt"] = np.tile(inst_data["MSL_alt"], (n_prof, 1))
                else:
                    store["Alt"] = inst_data["Alt"]
            elif name == "SourceProf":
                n_prof, n_lev = np.shape(store["Lat"])[:2]
                profile = np.arange(1, n_prof + 1, dtype=np.float32)  # 1-based
                store["SourceProf"] = np.tile(profile.reshape(-1, 1), (1, n_lev))
            elif name == "SourceFile":
                store["SourceFile"] = np.ones_like(store["SourceProf"]) * file_count
            else:
                store[name] = inst_data[name]
        return conform_vars(store, settings)

    return loop_days(settings, var_names, find_file, read_file, instrument, inst_info)
