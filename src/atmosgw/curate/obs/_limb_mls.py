"""
Limb-sounder reader for MLS temperature (NASA Level 2 HDF-EOS5 files,
``Temperature`` swath). Port of ``module_load_MLS.m`` (uobGravWaves/
MultiInstrumentFunctions, commit 78eef78), with the parts of ``get_MLS.m``
(corwin365/MatlabFunctions, commit c943399) that it uses.

Files: one per day, the file below the MLS folder whose path contains
``_<yyyy>d<ddd>``. If several match (e.g. v05-04 and v06-04 of the same day,
as for many days of 2026 in the group's MLS folder), the newest version is
used (``_limb_tools.newest_mls_version``), and the day is listed in the
output's ``multiple_matches`` attribute. The MATLAB loader does the same
since 2026-11-25; before, both took the first match in name order, there the
older version (Q8 in ``docs/quirks.md``).

The fields ``get_MLS`` reads are available: data fields
``L2gpValue`` (temperature), ``L2gpPrecision``, ``Quality``, ``Status``,
``Convergence``; geolocation fields ``ChunkNumber``, ``Latitude``,
``LineOfSightAngle``, ``LocalSolarTime``, ``Longitude``, ``OrbitGeodeticAngle``,
``Pressure``, ``SolarZenithAngle``, ``Time`` (seconds since 1993-01-01).

Each field keeps its stored number type (as ``get_MLS``), except that a field
containing its ``_FillValue`` is converted to float64 with those values NaN.
Negative ``L2gpPrecision`` values are meaningful (a priori influence) and are
kept. Heights are estimated from the pressure levels with
``atmosphere.pressure_to_height`` (U.S. Standard Atmosphere 1976), in double
precision (MLS pressures are float32; until 2026-11-25 the heights were computed
in single precision, as the old MATLAB ``p2h`` did).
"""

import h5py
import numpy as np

from atmosgw.curate.atmosphere import pressure_to_height
from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs._common import (
    _append_note,
    choose_file,
    day_label,
    find_in,
    note,
    note_sources,
    note_version,
    record,
    year_and_day,
)
from atmosgw.curate.obs._limb_tools import conform_vars, loop_days, newest_mls_version
from atmosgw.curate.readers import read_attrs

_SWATH = "Temperature"
DATA_FIELDS = ("L2gpValue", "L2gpPrecision", "Quality", "Status", "Convergence")
REQUIRED = ("L2gpValue", "Latitude", "Longitude", "Pressure", "Time")
GEOLOCATION_FIELDS = ("ChunkNumber", "Latitude", "LineOfSightAngle", "LocalSolarTime",
                      "Longitude", "OrbitGeodeticAngle", "Pressure", "SolarZenithAngle",
                      "Time")


def load(settings, inst_info, var_names):
    """Read all days; returns (data, file_list). See limb.py for the arguments."""
    note_sources(inst_info, {"Temp": "L2gpValue", "Lat": "Latitude",
                             "Lon": "Longitude", "Pres": "Pressure", "Time": "Time"})
    note(inst_info, "Alt", "estimated from Pressure with "
                           "atmosphere.pressure_to_height (U.S. Standard "
                           "Atmosphere 1976 [R1])")
    note(inst_info, "Time", "file Time: TAI seconds since 1993-01-01, converted "
                            "without leap seconds")

    def find_file(day):
        year, doy = year_and_day(day)
        found = find_in(inst_info, f"_{year:04d}d{doy:03d}")
        # several versions of a day: the newest (Q8; MATLAB takes the first)
        return choose_file(found, inst_info, day_label(day), newest_mls_version)

    def read_file(path, file_count, day):
        working = _get_mls(path)
        note_version(inst_info, path, pge_version(path))
        store = {}
        screened = None
        for name in var_names:
            if name in ("InstrumentBad", "QcIncomplete"):  # (data states: screen())
                if screened is None:
                    screened = screen(working, pge_version(path), inst_info)
                store[name] = screened[0 if name == "InstrumentBad" else 1]
            elif name == "Temp":
                store["Temp"] = working["L2gpValue"]
            elif name == "Lat":
                store["Lat"] = working["Latitude"]
            elif name == "Lon":
                store["Lon"] = working["Longitude"]
            elif name == "Alt":
                store["Alt"] = pressure_to_height(working["Pressure"])
            elif name == "Pres":
                store["Pres"] = working["Pressure"]
            elif name == "Time":
                store["Time"] = datenum(1993, 1, 1, 0, 0, working["Time"])
            elif name == "SourceProf":
                n = store["Lat"].size
                store["SourceProf"] = np.arange(1.0, n + 1)[:, None]  # 1-based
            elif name == "SourceFile":
                store["SourceFile"] = np.ones_like(store["SourceProf"]) * file_count
            elif name in working:
                store[name] = working[name]
            else:
                raise KeyError(f"Variable {name} not found in {path}")

        # per-profile and per-level variables onto [profiles x levels]
        n_prof, n_lev = store["Temp"].shape
        for name in ("Lat", "Lon", "Time", "SourceProf", "SourceFile"):
            if name not in store:  # (source indices only if asked for)
                continue
            store[name] = np.tile(store[name], (1, n_lev))
        store["Pres"] = np.tile(store["Pres"].T, (n_prof, 1))
        store["Alt"] = np.tile(store["Alt"].T, (n_prof, 1))
        for name in settings["additional_vars"]:
            field = store[name]
            if field.size == n_lev:
                store[name] = np.tile(field.reshape(1, -1, order="F"), (n_prof, 1))
            elif field.size == n_prof:
                store[name] = np.tile(field.reshape(-1, 1, order="F"), (1, n_lev))
        return conform_vars(store, settings)

    return loop_days(settings, var_names, find_file, read_file, "MLS", inst_info)


def _get_mls(path):
    """The fields get_MLS reads, in their stored types, 1-D fields as columns;
    a field containing its _FillValue becomes float64 with NaN there."""
    working = {}
    try:
        f = h5py.File(path, "r")
    except OSError as error:
        raise OSError(f"Cannot read {path}: {error}") from error
    with f:
        for group, names in (("Data Fields", DATA_FIELDS),
                             ("Geolocation Fields", GEOLOCATION_FIELDS)):
            for name in names:
                key = f"HDFEOS/SWATHS/{_SWATH}/{group}/{name}"
                if key not in f:
                    if name in REQUIRED:
                        raise KeyError(f"{key} not found in {path}")
                    working[name] = np.empty((0, 0))
                    continue
                dataset = f[key]
                value = dataset[()]
                fill = dataset.attrs.get("_FillValue")
                if fill is not None:
                    bad = value == np.asarray(fill, dtype=value.dtype).ravel()[0]
                    if bad.any():
                        value = value.astype(np.float64)
                        value[bad] = np.nan
                working[name] = value[:, None] if value.ndim == 1 else value
    return working


# The screening of the MLS data quality documents for temperature, the same
# in v5.0x (Livesey et al., MLS Version 5.0x Level 2 and 3 data quality and
# description document, section 3.22.6) and v6.x (section 3.23.6), checked
# against both documents on 2026-11-26. Not applied: cloud screening with
# the ice water content product (a separate swath), and the optional
# precision tests for 261-178 hPa.
USEFUL_RANGE_HPA = (0.00046, 261.0)  # "not recommended" outside
QUALITY_LOW_P = 0.2  # at 83 hPa and smaller pressures: Quality must exceed
QUALITY_HIGH_P = 0.9  # at 100 hPa and larger pressures
QUALITY_SPLIT_HPA = 90.0  # (between the 100 and 83 hPa levels)
CONVERGENCE_MAX = 1.03  # Convergence must be below
SCREENED_VERSIONS = ("v5.", "v6.")


def screen(working, version_text, inst_info):
    """
    ([profiles x levels] bad, [profiles x levels] incomplete), as floats 1.0
    or 0.0. bad: where the MLS data quality document says not to use a
    temperature value (see above). Profile rules (Status odd, Convergence >=
    1.03, Quality at or below the threshold for the level's pressure) mark
    the whole profile's affected levels; level rules (precision <= 0,
    pressure outside 261-0.00046 hPa) mark points. A rule whose field is
    missing from the file is not applied, the others are, and every value of
    the file is incomplete (``qc_incomplete``).
    Versions the documents do not cover are not screened (noted).
    """
    temp = np.asarray(working["L2gpValue"])
    shape = temp.shape
    none = np.zeros(shape)
    version = version_text.split()[1] if version_text else ""
    if not version.lower().startswith(SCREENED_VERSIONS):
        note(inst_info, "Temp", f"no MLS screening rules for version "
                                f"{version or 'unknown'}: not screened")
        return none, np.zeros(shape)
    missing = [f for f in ("Status", "Quality", "Convergence", "L2gpPrecision")
               if np.size(working[f]) == 0]
    note(inst_info, "Temp", "MLS data quality document screening applied "
                            "(v5.0x section 3.22.6, v6.x section 3.23.6): "
                            "Status even, Convergence < 1.03, Quality > 0.2 at "
                            "83 hPa and above / > 0.9 at 100 hPa and below, "
                            "precision > 0, 261-0.00046 hPa; failures are "
                            "instrument_flag (no cloud screening)")
    if missing:  # (appended: the screening note above holds the Temp key)
        text = (f"MLS quality fields missing from a file ({', '.join(missing)}): "
                "their rules not applied there, its values qc_incomplete")
        notes = record(inst_info)["notes"]
        if text not in notes.get("Temp", ""):
            _append_note(notes, "Temp", text)
    pres = np.asarray(working["Pressure"], dtype=float).ravel()[None, :]
    with np.errstate(invalid="ignore"):
        # (the files' levels are float32 and nominal: 261 hPa is stored as
        # 261.016, so the range is compared to within 1e-4)
        bad = ((pres < USEFUL_RANGE_HPA[0] * (1 - 1e-4))
               | (pres > USEFUL_RANGE_HPA[1] * (1 + 1e-4)))
        if "Status" not in missing:
            status = np.asarray(working["Status"], dtype=float).reshape(-1, 1)
            bad = bad | ~np.isfinite(status) | (np.mod(status, 2) == 1)
        if "Convergence" not in missing:
            convergence = np.asarray(working["Convergence"], dtype=float).reshape(-1, 1)
            bad = bad | ~(convergence < CONVERGENCE_MAX)
        if "Quality" not in missing:
            quality = np.asarray(working["Quality"], dtype=float).reshape(-1, 1)
            threshold = np.where(pres < QUALITY_SPLIT_HPA, QUALITY_LOW_P,
                                 QUALITY_HIGH_P)
            bad = bad | ~(quality > threshold)
        if "L2gpPrecision" not in missing:
            bad = bad | ~(np.asarray(working["L2gpPrecision"], dtype=float) > 0)
    incomplete = np.ones(shape) if missing else none
    return np.broadcast_to(bad, shape).astype(float), incomplete


def pge_version(path):
    """'PGEVersion <v> (DOI <doi>)' from an Aura HDF-EOS5 file's attributes."""
    attrs = read_attrs(path, "HDFEOS/ADDITIONAL/FILE_ATTRIBUTES")
    version = attrs.get("PGEVersion", "")
    doi = attrs.get("identifier_product_doi") or attrs.get("identifier_product_DOI")
    if not version:
        return ""
    return f"PGEVersion {version}" + (f" (DOI {doi})" if doi else "")
