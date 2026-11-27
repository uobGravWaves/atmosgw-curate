"""
Limb-sounder reader for output of the Wright and Hindley (2018) sampler
(https://github.com/corwin365/Model-Sampling).
Port of ``module_load_sampler.m`` (uobGravWaves/MultiInstrumentFunctions,
commit 78eef78).

File format: one MATLAB .mat file per day, ``<sample_path>/sampled_<datenum>.mat``
(e.g. ``sampled_733043.mat`` for 2007-01-01), holding a struct
``Sampled_Data`` whose fields are [levels x profiles] (transposed on reading):
``Alt``, ``Lat``, ``Lon``, ``Time``, ``T`` (temperature), ``Prs`` (pressure)
and ``SourceProf``.

Where a day has no such file, the first netCDF file below ``sample_path``
whose name ends ``<yyyy>d<ddd>.nc`` is read instead (new, not in the MATLAB):
the sampler's netCDF output as in the SWANS DYAMOND runs, e.g.
``HIRDLS-dyamond_geos3km-2020d020.nc``, holding the same variables (``Temp``
for temperature, ``ObsProf`` for ``SourceProf``) as [levels x profiles],
-9999 (``_FillValue``) as NaN. ``Time`` is a datenum in both.

Change log:
2026/11/27 [Claude] netCDF sampler output (SWANS) read where there is no .mat file
"""

import os

import numpy as np

from atmosgw.curate.obs._common import (
    choose_file,
    day_label,
    find_in,
    note_sources,
    year_and_day,
)
from atmosgw.curate.obs._limb_tools import conform_vars, loop_days
from atmosgw.curate.readers import read_mat, read_netcdf

# netCDF sampler output (SWANS): output names -> file names
NC_NAMES = {"Temp": "Temp", "Pres": "Prs", "Lat": "Lat", "Lon": "Lon", "Alt": "Alt",
            "Time": "Time", "SourceProf": "ObsProf"}


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    sample_path = settings["sample_path"]
    if sample_path is None or not os.path.isdir(sample_path):
        raise FileNotFoundError(f"{sample_path} does not exist")
    note_sources(inst_info, {"Temp": "Sampled_Data.T", "Pres": "Sampled_Data.Prs",
                             "Lat": "Sampled_Data.Lat", "Lon": "Sampled_Data.Lon",
                             "Alt": "Sampled_Data.Alt", "Time": "Sampled_Data.Time"})

    inst_info["path"] = str(sample_path)  # (searched for netCDF output)

    def find_file(day):
        path = os.path.join(str(sample_path), f"sampled_{day}.mat")
        if os.path.isfile(path):
            return path
        year, doy = year_and_day(day)
        found = [f for f in find_in(inst_info, f"*{year:04d}d{doy:03d}.nc")
                 if f.endswith(f"{year:04d}d{doy:03d}.nc")]
        return choose_file(found, inst_info, day_label(day)) or path

    def read_file(path, file_count, day):
        if path.endswith(".nc"):
            names = [NC_NAMES.get(n, n) for n in var_names if n != "SourceFile"]
            try:
                working = read_netcdf(path, names)
            except KeyError as error:
                raise KeyError(f"Variable {error} not found in {path}") from error
            inst_data = {n: working[NC_NAMES.get(n, n)] for n in var_names
                         if n != "SourceFile"}
            note_sources(inst_info, {n: NC_NAMES[n] for n in NC_NAMES
                                     if n != "SourceProf"})
        else:
            inst_data = dict(read_mat(path)["Sampled_Data"])
            inst_data["Temp"] = inst_data.pop("T")
            inst_data["Pres"] = inst_data.pop("Prs")
        store = {}
        for name in var_names:
            if name == "SourceFile":
                store["SourceFile"] = np.ones_like(store["SourceProf"]) * file_count
            elif name not in inst_data:
                raise KeyError(f"Variable {name} not found in {path}")
            else:
                store[name] = np.asarray(inst_data[name]).T
        return conform_vars(store, settings)

    return loop_days(settings, var_names, find_file, read_file, "Sampler", inst_info)
