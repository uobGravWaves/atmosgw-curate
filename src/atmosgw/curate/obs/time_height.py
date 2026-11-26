"""
Load ground-based time-height data (lidars, later radars) in a common format.

Load the measurements of an instrument at a fixed site as one [time x level]
grid, in the common "time_height" format (``docs/design/data_formats.md``),
so that time series of profiles are not treated as satellite profiles.
First instrument: the ALOMAR lidar (Andoya, Norway). The format also defines
wind variables, for radars, but no radar reader is written yet.

``load_limb_sounders(..., "ALOMAR")`` still works (the MATLAB-parity path,
port of ``module_load_ALOMAR.m``); it warns and points here.

Inputs and outputs
------------------
``load_time_height(time_range, instrument, **options)`` -> ``xarray.Dataset``

``time_range``: one time (that day) or two (start and end), as
``numpy.datetime64``, ISO strings or ``datetime`` objects, as for
``load_limb_sounders`` (numbers are refused).

``instrument``: ``"ALOMAR"`` (any case).

Options (type, default), with the same meaning as in ``load_limb_sounders``:

- ``time_handling`` (str, "fuzzy"): "strict", "generous" or "fuzzy". Each
  row (time step) is selected by its own time.
- ``height_range`` ((min, max), km, (0, inf)): levels outside it are dropped.
- ``keep_outliers`` (bool, False): keep values outside the physical range
  (their state is still set).
- ``additional_vars`` (list of str, ()): extra variables, by their names in
  the files (ALOMAR: ``T2_lidar``, ``T2_lidar_unc``, ``T_ecmwf``). [times x
  heights] arrays become (time, level), [times] arrays (time). They are
  passed through as read: outliers in temperature are not removed from them.
- ``file_source`` (bool, False): also return ``source_profile``,
  ``source_file`` and ``original_files``.
- ``file_path`` (str or Path, None): read the files from this folder.
- ``verbose`` (bool, False), ``date_warning`` (bool, True).

The other options of ``load_limb_sounders`` (``original_z``,
``height_scale``, ``interp_type``, ``extrapolate``, ``lat_range``,
``lon_range``, ``hindley_pws``, ``ignore_non_modal``, ``dx``, ``misc_path``,
``misc_identifier``, ``sample_path``) do not apply and raise TypeError: the
data are always on their own heights, at one site.

**Output**: an ``xarray.Dataset`` with dimensions ``time`` and ``level``
(and ``file`` with ``file_source``):

- ``time`` (time): ``datetime64[us]``, sorted (stable: rows with the same
  time stay in file order). Times can repeat (see Science);
- ``altitude`` (level): km, the union of the heights of the sessions in the
  output;
- ``latitude``, ``longitude`` (no dimensions): the site;
- ``temperature`` (time, level): K, with ``temperature_flag`` and
  ``temperature_uncertainty`` (as supplied by the files);
- ``session`` (time): 1-based measurement session (HDF5 group) in the
  output, in file order;
- ``source_profile`` (time): 1-based row through the yearly file, counting
  every group in name order (as the limb loader's); ``source_file`` (time):
  1-based index into ``original_files`` (``file_source``).

Data variables are float64, as in the files. States set (shared rules,
``data_formats.md`` section 6.3): ``no_observation`` (heights a session does
not have, and above or below a row's values), ``missing`` (NaN between a
row's first and last value), ``outside_physical_range`` (temperature < 100
K, or > 400 K below 100 km; removed unless ``keep_outliers``, together with
its uncertainty). ``availability`` lists ``pressure`` and the winds as
``not_in_files``. Attributes: as the other loaders (``format``
"time_height", ``format_version``, ``instrument``, ``product``,
``characterisation``, ``data``, ``preprocessing``, ``heights`` "original",
``site_name``, ``time_range``, ``availability``, ``source``,
``atmosgw_curate_version``, ``options``, ``files_read``,
``files_not_found``, ``multiple_matches``, ``state_counts``).

``data.to_profiles`` turns the output into the profile format (one profile
per row, pressure estimated from height).

Science
-------
- ALOMAR: ``temperature`` is ``T1_lidar`` (the higher-resolution product, as
  the limb loader), ``temperature_uncertainty`` ``T1_lidar_unc``. What kind
  of uncertainty it is (e.g. 1-sigma random) is not stated in the files
  (UNVERIFIED). Site 69.38 N, 16.08 E (as ``module_load_ALOMAR.m``).
- Heights: every session's heights lie on a 0.15 km lattice, so the output
  holds every value at its own height, without interpolation; heights are
  the lattice values, within 1e-11 km of those in the files.
- Repeated times: about 29% of ALOMAR rows (77,579 of 270,317 in the files
  of 2003-2026) share their time stamp with one to three other rows, with
  different values, mostly at overlapping heights. Why is not stated in the
  files (UNVERIFIED; two beams?). All rows are kept, as by the limb loader.
- Time stamps: what a stamp means (start or centre of the integration) is
  not stated in the files (UNVERIFIED).
- Memory: a whole archive (24 years) is about 1.7 GB per float64 variable.

Origin
------
New code. The ALOMAR file layout follows the port of ``module_load_ALOMAR.m``
(``_limb_alomar.py``).

Validation
----------
``tests/test_time_height_known.py`` (synthetic files: values, states,
options), ``tests/test_time_height_alomar.py`` (real files, public
Datasets: values, rows and states identical to ``load_limb_sounders(...,
"ALOMAR", original_z=True)`` on five date ranges; a wider run on 34 dates,
2003-2026, made once on 2026-11-26), ``tests/test_states_conformance.py``
and ``tests/test_data_convert.py``. Full record: ``docs/validation.md`` in
the repository.

Status
------
Proposal (format version 1 is a draft for agreement). AI-originated
(Claude), from the plan agreed with C. Wright on 2026-11-26; to be reviewed
by C. Wright.

Usage
-----
>>> from atmosgw.curate.obs import load_time_height  # doctest: +SKIP
>>> ds = load_time_height("2007-01-23", "ALOMAR")  # doctest: +SKIP
>>> ds.temperature.sel(time=slice("2007-01-23T13", "2007-01-23T14"))  # doctest: +SKIP

Change log
----------
2026/11/26 [Claude] initial version (ALOMAR)
2026/11/26 [Claude] docstring: Validation scope (dates, public Datasets)

"""

import os
import warnings

import numpy as np
import xarray as xr

from atmosgw.curate.data import states
from atmosgw.curate.data.formats import FORMATS
from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs import _limb_alomar, _th_alomar
from atmosgw.curate.obs._common import (
    COORD_ATTRS,
    check_additional_vars,
    check_flags,
    check_time_handling,
    expand_time_range,
    file_versions,
    find_instrument,
    load_attrs,
    outside_physical_range,
    record,
    to_datenums,
    warn_if_no_data,
    warn_missing_files,
    warn_mixed_versions,
)
from atmosgw.curate.obs.limb import _state_flags
from atmosgw.curate.paths import dataset_dir

_FAR = datenum(9999, 999, 999)  # "still running"

# instruments: dates with data, reader, site, product, and the format's
# variables that the files do not have
INSTRUMENTS = {
    "ALOMAR": {"time_range": (datenum(2003, 1, 1), _FAR), "reader": _th_alomar,
               "site_name": "ALOMAR, Andoya, Norway",
               "latitude": _limb_alomar.SITE_LAT, "longitude": _limb_alomar.SITE_LON,
               "product": "temperature",
               "not_in_files": ("pressure", "eastward_wind", "northward_wind",
                                "upward_wind")},
}

# options of load_limb_sounders that have no meaning here
_LIMB_ONLY = ("original_z", "height_scale", "interp_type", "extrapolate",
              "lat_range", "lon_range", "hindley_pws", "ignore_non_modal", "dx",
              "misc_path", "misc_identifier", "sample_path")

_UNCERTAINTY_COMMENT = ("as supplied in the files; what kind of uncertainty it "
                        "is (e.g. 1-sigma random) is not stated there (UNVERIFIED)")


def load_time_height(
    time_range,
    instrument,
    *,
    time_handling="fuzzy",
    height_range=(0, np.inf),
    keep_outliers=False,
    additional_vars=(),
    file_source=False,
    file_path=None,
    verbose=False,
    date_warning=True,
    **other,
):
    # (the module docstring, so that help(load_time_height) shows it all)
    if other:
        limb = [name for name in other if name in _LIMB_ONLY]
        if limb:
            raise TypeError(f"load_time_height does not take {', '.join(limb)} "
                            "(options of load_limb_sounders: time-height data are "
                            "always on their own heights, at one site).")
        raise TypeError("load_time_height got unexpected keyword argument(s): "
                        + ", ".join(other))
    options = {k: v for k, v in locals().items()
               if k not in ("time_range", "instrument", "other")}
    call = {"time_range": time_range, "instrument": instrument, **options}
    settings, inst_info, reader = _prepare(time_range, instrument, options)
    rows, file_list = reader.load(settings, inst_info)
    settings["record"] = record(inst_info)
    settings["files_read"] = list(file_list)
    warn_mixed_versions(settings, stacklevel=2)
    warn_missing_files(settings, stacklevel=2, what="yearly files")
    if (settings["scalar_time"] and settings["time_handling"] == "strict"
            and rows["time_posix"].size == 0):
        warnings.warn(
            "A single time_range value with time_handling='strict' selects only "
            "times exactly at that instant, and none were found. Use "
            "time_handling 'generous' or 'fuzzy' (default) to get the whole day.",
            stacklevel=2,
        )
    return _to_dataset(rows, file_list, settings, call)


# supporting functions


def _prepare(time_range, instrument, options):
    """Check the arguments; return (settings, inst_info, reader)."""
    check_flags(options, ("keep_outliers", "file_source", "verbose", "date_warning"))
    check_time_handling(options)
    check_additional_vars(options, {
        "temperature", "temperature_uncertainty", "session", "time", "altitude",
        "latitude", "longitude", "source_profile", "source_file", "original_files",
        "original_file_version"})
    name = find_instrument(instrument, INSTRUMENTS,
                           " (satellite profiles: load_limb_sounders)")
    info = INSTRUMENTS[name]

    height_range = tuple(float(h) for h in options["height_range"])
    if len(height_range) != 2 or not all(h >= 0 for h in height_range):
        raise ValueError("height_range must be two heights >= 0 (km).")

    raw = to_datenums(time_range)
    if options["date_warning"]:
        warn_if_no_data(raw, name, info["time_range"], stacklevel=3)
    expanded, scalar = expand_time_range(raw, options["time_handling"])
    settings = dict(options, instrument=name, time_range=expanded,
                    scalar_time=scalar, height_range=height_range,
                    additional_vars=list(options["additional_vars"]))

    inst_info = {"path": "", "inst": name}
    if options["file_path"] is not None:
        inst_info["path"] = str(options["file_path"])
    else:
        inst_info["path"] = str(dataset_dir(name, "raw", must_exist=False))
    if not os.path.isdir(inst_info["path"]):
        given = ("file_path" if options["file_path"] is not None
                 else "the data root and datasets.toml")
        raise FileNotFoundError(f"The {name} data folder {inst_info['path']} does not "
                                f"exist (from {given}; see atmosgw.curate.paths).")
    settings["data_path"] = inst_info["path"]
    return settings, inst_info, info["reader"]


def _to_dataset(rows, file_list, settings, call):
    """The reader's rows as a Dataset in the time_height format."""
    name = settings["instrument"]
    info = INSTRUMENTS[name]
    temp_name = f"{_limb_alomar.RESOLUTION}_lidar"
    unc_name = f"{_limb_alomar.RESOLUTION}_lidar_unc"
    rec = settings["record"]

    # heights: lattice values, then only those in height_range
    altitude = np.round(rows["level_index"] * _th_alomar.LEVEL_STEP, 2)
    low, high = min(settings["height_range"]), max(settings["height_range"])
    levels = np.flatnonzero((altitude >= low) & (altitude <= high))
    altitude = altitude[levels]

    # rows in time order (stable: repeated times keep file order)
    order = np.argsort(rows["time_posix"], kind="stable")

    def per_point(key):
        return rows[key][order][:, levels]

    temp, unc = per_point(temp_name), per_point(unc_name)
    present = per_point("present")
    masks = {"outside_physical_range": outside_physical_range(temp, altitude[None, :]),
             "no_observation": ~present}
    if not settings["keep_outliers"]:
        temp[masks["outside_physical_range"]] = np.nan
        unc[masks["outside_physical_range"]] = np.nan
    flag = _state_flags(temp, temp.shape, masks)
    unc[np.isnan(temp)] = np.nan  # (no flag of its own: NaN where temperature is)

    # POSIX seconds to the microsecond
    micro = np.round(rows["time_posix"][order] * 1e6).astype(np.int64)
    time = micro.astype("datetime64[us]")

    tl = ("time", "level")
    variables = {
        "temperature": (tl, temp, {
            "units": "K", "standard_name": "air_temperature",
            "source_variable": rec["source_variables"]["Temp"],
            "origin": "retrieved",
            "ancillary_variables": "temperature_flag temperature_uncertainty"}),
        "temperature_flag": (tl, flag, states.flag_attrs()),
        "temperature_uncertainty": (tl, unc, {
            "units": "K", "long_name": "temperature uncertainty",
            "source_variable": rec["source_variables"]["TempUnc"],
            "comment": _UNCERTAINTY_COMMENT}),
        "session": (("time",), rows["session"][order].astype(np.int32), {
            "long_name": "measurement session (HDF5 group) in this output, "
                         "1-based, in file order"}),
    }
    for extra in settings["additional_vars"]:
        value = rows[extra][order]
        dims = ("time",) if value.ndim == 1 else tl
        if value.ndim == 2:
            value = value[:, levels]
        variables[extra] = (dims, value, {"source_variable": extra})
    coords = {
        "time": (("time",), time, {
            **COORD_ATTRS["time"],
            "source_variable": rec["source_variables"]["Time"],
            "comment": "POSIX seconds in the files, kept to the microsecond; "
                       "start or centre of the integration not stated "
                       "(UNVERIFIED); some times repeat"}),
        "altitude": (("level",), altitude, {
            **COORD_ATTRS["altitude"],
            "source_variable": rec["source_variables"]["Alt"],
            "comment": f"union of the sessions' heights, on their common "
                       f"{_th_alomar.LEVEL_STEP} km grid (within 1e-11 km of "
                       "the file values)"}),
        "latitude": ((), float(info["latitude"]), {
            **COORD_ATTRS["latitude"], "comment": rec["notes"]["Lat"]}),
        "longitude": ((), float(info["longitude"]), {
            **COORD_ATTRS["longitude"], "comment": rec["notes"]["Lon"]}),
    }
    if settings["file_source"]:
        variables["source_profile"] = (("time",), rows["source_profile"][order]
                                       .astype(np.int32), {
            "long_name": "row number in the source file (1-based, counting "
                         "every session in name order)"})
        variables["source_file"] = (("time",), rows["source_file"][order]
                                    .astype(np.int32), {
            "long_name": "index into original_files (1-based)"})
        coords["original_files"] = (("file",), np.array(file_list, dtype=str))
        coords["original_file_version"] = (("file",), np.array(
            file_versions(rec, file_list), dtype=str),
            {"long_name": "data version stated by each file"})

    spec = FORMATS["time_height"]["variables"]
    availability = {n: ("not_in_files" if n in info["not_in_files"]
                        else "not_requested")
                    for n, s in spec.items() if s.get("tracked")
                    and n not in variables and n not in coords}
    common = load_attrs(settings, call, "atmosgw.curate.obs.load_time_height")
    product = info["product"]
    attrs = {
        "format": "time_height",
        "format_version": FORMATS["time_height"]["version"],
        "instrument": name, "product": product,
        "characterisation": f"{name}-{product.upper()}-placeholder/0",
        "data": "raw", "preprocessing": "", "heights": "original",
        "site_name": info["site_name"],
        "time_range": common.pop("time_range"),
        "availability": states.format_availability(availability),
        **common,
        "state_counts": states.counts_text({"temperature": flag}),
    }
    return xr.Dataset(variables, coords=coords, attrs=attrs)


load_time_height.__doc__ = __doc__

__all__ = ["load_time_height", "INSTRUMENTS"]
