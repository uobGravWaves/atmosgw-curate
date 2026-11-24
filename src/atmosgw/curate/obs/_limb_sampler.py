"""
Limb-sounder reader for output of the Wright and Hindley (2018) sampler.
Port of ``module_load_sampler.m`` (uobGravWaves/MultiInstrumentFunctions,
commit 78eef78).

File format: one MATLAB .mat file per day, ``<sample_path>/sampled_<datenum>.mat``
(e.g. ``sampled_733043.mat`` for 2007-01-01), holding a struct
``Sampled_Data`` whose fields are [levels x profiles] (transposed on reading):
``Alt``, ``Lat``, ``Lon``, ``Time``, ``T`` (temperature), ``Prs`` (pressure)
and ``SourceProf``.
"""

import os

import numpy as np

from atmosgw.curate.obs._limb_tools import conform_vars, loop_days
from atmosgw.curate.readers import read_mat


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    sample_path = settings["sample_path"]
    if sample_path is None or not os.path.isdir(sample_path):
        raise FileNotFoundError(f"{sample_path} does not exist")

    def find_file(day):
        return os.path.join(str(sample_path), f"sampled_{day}.mat")

    def read_file(path, file_count, day):
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

    return loop_days(settings, var_names, find_file, read_file, "Sampler")
