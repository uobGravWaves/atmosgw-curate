"""
Nadir-sounder reader for AIRS 2-D brightness temperatures (L. Hoffmann's
granule files, with perturbations and variances).

Files, below the AIRS 2-D folder (``datasets.toml``, ``AIRS`` ``2d``):
``YYYY/DDD/airs_YYYY_DDD_GGG.nc``, or ``DDD/airs_YYYY_DDD_GGG.nc`` (copies
holding one year), else that name anywhere below the folder, one per
granule: 135 along-track x 90 cross-track points (dimensions ``NTRACK``,
``NXTRACK``), with ``time`` (seconds since 2000-01-01T00:00Z), ``lat``,
``lon`` and, for each channel group, ``bt_<group>`` (K) and, except 8.1
micron, its perturbation ``bt_<group>_pt`` (K) and variance
``bt_<group>_var`` (K^2).

Channel groups (``channels`` option; default all, in this order):

========== ==================== ============================== ============
name       file variable        what (from the files)           nominal
                                                                cm-1
========== ==================== ============================== ============
4mu        bt_4mu               4.3 micron                      2325.6
8mu        bt_8mu               8.1 micron (no perturbation)    1234.6
15mu_low   bt_15mu_low          15 micron, low altitudes         666.7
15mu_high  bt_15mu_high         15 micron, high altitudes        666.7
========== ==================== ============================== ============

The wavenumbers are nominal (10^4 / wavelength in micron, from the files'
descriptions), not the channels actually averaged by the producer:
UNVERIFIED. Each channel is labelled in ``channel_name``.

Variables (reader contract, ``docs/design/nadir_loader.md`` section 5):
``Lat``, ``Lon``, ``Time`` (datenums), ``Temp`` (brightness temperature)
[cross x along x channel], ``Channel`` (wavenumbers), ``ChannelName``, per row
``Granule``, ``Day``, ``SourceFile``, ``SourceRow``, and the producer's
perturbation and variance ``BTPert``, ``BTVar`` [cross x along x channel]
(NaN for 8.1 micron), declared in ``DerivedVars`` so the loader regrids them
and returns them as ``brightness_temperature_perturbation`` and
``brightness_temperature_variance`` with ``origin`` "derived".

Granules are chosen as for the 3-D AIRS reader (``_nadir_airs``): by the
mean of their times, or given explicitly; missing granules are recorded.
"""

import os

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
    record,
    year_and_day,
)
from atmosgw.curate.obs._nadir_airs import (
    FIRST_START,
    GRANULE_DAYS,
    J2000,
    candidate_granules,
)
from atmosgw.curate.readers import read_netcdf

# name: (file variable, nominal wavenumber cm-1, description)
CHANNELS = {
    "4mu": ("bt_4mu", 1e4 / 4.3, "4.3 micron"),
    "8mu": ("bt_8mu", 1e4 / 8.1, "8.1 micron"),
    "15mu_low": ("bt_15mu_low", 1e4 / 15, "15 micron, low altitudes"),
    "15mu_high": ("bt_15mu_high", 1e4 / 15, "15 micron, high altitudes"),
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


def find_file(root, day, granule, inst_info=None):
    """Path of a granule file (YYYY/DDD/ or DDD/ below root, else, with
    inst_info, anywhere below root), or None."""
    year, doy = year_and_day(day)
    name = f"airs_{year:04d}_{doy:03d}_{granule:03d}.nc"
    for folder in (os.path.join(root, f"{year:04d}", f"{doy:03d}"),
                   os.path.join(root, f"{doy:03d}")):
        path = os.path.join(folder, name)
        if os.path.isfile(path):
            return path
    if inst_info is not None:
        return locate(inst_info, os.path.join(root, name), name, day_label(day), root)
    return None


def load(settings, inst_info, var_names):
    """Read the granules; returns (data, file_list). See the module docstring."""
    channels = list(settings.get("channels") or CHANNELS)
    unknown = [c for c in channels if c not in CHANNELS]
    if unknown:
        raise ValueError(f"Unknown AIRS 2-D channel(s) {unknown}; choose from "
                         f"{', '.join(CHANNELS)}.")
    note_sources(inst_info, {"Lat": "lat", "Lon": "lon", "Time": "time",
                             "Temp": ", ".join(CHANNELS[c][0] for c in channels)})
    note(inst_info, "Time", "time: seconds since 2000-01-01T00:00Z")
    explicit = settings.get("granules") is not None
    t0, t1 = (-np.inf, np.inf) if explicit else settings["time_range"]
    rec = record(inst_info)

    granules, file_list, fill_counts = [], [], {}
    for day, g in candidate_granules(settings):
        path = find_file(inst_info["path"], day, g, inst_info)
        if path is None:
            mid = day + FIRST_START + (g - 0.5) * GRANULE_DAYS
            r0, r1 = (-np.inf, np.inf) if explicit else settings.get(
                "request_range", settings["time_range"])
            if r0 <= mid <= r1:  # (as the 3-D reader)
                rec["files_not_found"].append(f"{day_label(day)} granule {g}")
            continue
        names = ["time", "lat", "lon"]
        for c in channels:
            base = CHANNELS[c][0]
            names += [base] + [base + d["suffix"] for d in DERIVED.values()]
        values = _read_present(path, names)
        # fill values and bad times, as the 3-D readers (_common.nadir_fill_to_nan)
        span = ((day - J2000) * 86400, (day + 1 - J2000) * 86400)
        counts = nadir_fill_to_nan(values["lat"], values["lon"], values["time"],
                                   np.zeros(values["lat"].shape + (1,)), span)
        for name in names[3:]:
            if name in values:
                fill = np.isin(values[name], FILL_VALUES)
                counts["temp"] += int(fill.sum())
                values[name][fill] = np.nan
        for key, n in counts.items():
            fill_counts[key] = fill_counts.get(key, 0) + n
        time = J2000 + values["time"].T / 86400  # [cross x along]
        mean = np.nanmean(time) if np.any(np.isfinite(time)) else np.nan
        if not explicit and not (t0 <= mean <= t1):
            continue
        file_list.append(os.path.basename(path))
        n_along = values["lat"].shape[0]

        def stack(suffix, values=values, n_along=n_along):
            planes = []
            for c in channels:
                v = values.get(CHANNELS[c][0] + suffix)
                planes.append(np.full(values["lat"].T.shape, np.nan) if v is None
                              else v.T)
            return np.stack(planes, axis=-1)

        part = {"Lat": values["lat"].T, "Lon": values["lon"].T, "Time": time,
                "Temp": stack(""), "Granule": np.full(n_along, float(g)),
                "Day": np.full(n_along, float(day)),
                "SourceFile": np.full(n_along, float(len(file_list))),
                "SourceRow": np.arange(1, n_along + 1, dtype=float)}
        for key, d in DERIVED.items():
            part[key] = stack(d["suffix"])
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
        data = {name: np.empty((90, 0, len(channels))) if name in ("Temp", *DERIVED)
                else (np.empty((90, 0)) if name in ("Lat", "Lon", "Time")
                      else np.empty(0)) for name in names}
    data["Channel"] = np.array([CHANNELS[c][1] for c in channels])
    data["ChannelName"] = np.array(channels)
    data["DerivedVars"] = DERIVED
    if "8mu" in channels:
        note(inst_info, "BTPert", "no perturbation or variance in the files for 8mu "
                                  "(NaN)")
    return data, file_list


def _read_present(path, names):
    """The variables of names that the file has (others left out)."""
    with netCDF4.Dataset(path) as f:
        present = [n for n in names if n in f.variables]
    return read_netcdf(path, present)
