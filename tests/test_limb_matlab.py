"""Regression tests: load_limb_sounders against the MATLAB get_limbsounders.

Reference: <reference_dir>/get_limbsounders_reference.mat, made 2026-11-25 by
corwin365/MatlabFunctions-Tests@1c89c4b (test_get_limbsounders.m; 37 cases)
with MATLAB R2024b Update 6, from uobGravWaves/MultiInstrumentFunctions@52518fc
and corwin365/MatlabFunctions@3b84a6b. Its 33 earlier cases are bit-identical
to the previous reference, kept beside it as
get_limbsounders_reference_9db1b4b.mat (made from
MultiInstrumentFunctions@9db1b4b and MatlabFunctions@d54b177); older ones:
get_limbsounders_reference_9eb6bde.mat (before Q1-Q3) and
get_limbsounders_reference_78eef78.mat (before Q6). Its SHA-256 is
pinned below; a different file fails the test rather than being used.

Synthetic inputs: <reference_dir>/get_limbsounders_synthetic, written by
tests/matlab/make_limb_synthetic_inputs.m (checked 2026-11-24 to reproduce all
14 synthetic reference cases in MATLAB; rewritten 2026-11-25 with the new
h2p.m, which changed only the sampler file's pressure; the previous files are
kept in get_limbsounders_synthetic_20261124).

What is compared: the private layer (_load_limb_arrays), whose output has the
MATLAB names and datenum times, field by field: the same fields, shapes and NaN
pattern, and equal values.

Tolerances: bit-identical, except:
- heights computed from pressure (MLS, heights kept): 1e-15 relative (log()
  rounding, as pressure below); largest seen 2.0e-16 (2026-11-25 sweep);
- pressure computed from height (Misc and PW files have no pressure): 1e-15
  relative, because exp() and log() round differently in the last bit (see
  tests/test_atmosphere.py); largest seen 4.9e-16; above 86 km, where the
  standard atmosphere's table is interpolated in log(pressure), 1e-14
  relative (ALOMAR_origz_high; largest seen 2.0e-15);
- pressure when interpolating to the height scale: 1e-14 relative, because it
  is interpolated in log(pressure) (Q6 fix) and log() and exp() round
  differently in the last bit; largest seen 2.0e-15 (2026-11-25);
- HIRDLS latitude and longitude: 1e-12 degrees absolute, because the
  geolocation correction uses azimuth and reckon, whose trig functions round
  differently in the last bits (see tests/test_geodesy.py); largest seen
  7.1e-14 degrees (2026-11-24);
- dx where positions differ as above: 1e-12 relative; largest seen 2e-14;
  dx of other instruments (haversine trig): the same; largest seen 1.2e-13
  (2026-11-25 sweep);
- ALOMAR times (POSIX seconds converted through datetime in MATLAB): 1.2e-10
  days absolute, one unit in the last place of a datenum (10 microseconds).
NaN patterns, shapes and the number of profiles must match exactly.

Q12 (boxes need one level inside, both coordinates): the MATLAB rows of a
boxed case are reduced to those (_q12_rows) before the comparison.

Deliberate differences from MATLAB since 0.5.1 (docs/quirks.md; the last
release that reproduces MATLAB is 0.5.1): the MATLAB output of the affected
cases is changed as documented before the comparison (_python_changes):
ACE pressure x 1013.25 (Q9); GNSS pressure / 100 (Q9), positions from
MetaData for profiles with none (Q10; then 1e-12 degrees, as the longitude
goes through unwrap and wrap), times from the occultation id where the file's
is hour-rounded (Q10), and only the rows Python keeps (Q11: one copy per
occultation; a profile given a MetaData position outside a box is left out).
The rows Python leaves out are checked to be exactly those.

Real-data cases use the MATLAB test's dates: 2007-01-01 (and -02) for most;
2007-01-02 for MIPAS, 2007-06-01 for SOFIE and 2007-01-23 for ALOMAR; and,
added after the sweep of tests/matlab/sweep (2026-11-25), 2005-01-15 for
ALOMAR with heights to 120 km and 2026-01-01 for MLS (v05-04 and v06-04 files
on the same day).

Skipped (not passed) if reference_dir is not set or the files are absent.
"""

import hashlib
import warnings
from pathlib import Path

import numpy as np
import pytest

from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs.limb import _load_limb_arrays
from atmosgw.curate.paths import DataDirNotSetError, dataset_dir, reference_dir
from atmosgw.curate.readers import read_mat

REFERENCE_SHA256 = "7a551d0650b5d9d6b0869905ec1fbf31bdd12945789e03f2c35935c8aefe08a7"
PRESSURE_RTOL = 1e-15
POSITION_ATOL = 1e-12  # degrees
DX_RTOL = 1e-12
LOG_PRESSURE_RTOL = 1e-14  # interpolated pressure (log/exp)
HEIGHT_RTOL = 1e-15  # heights computed from pressure (MLS, heights kept)
HIGH_PRESSURE_RTOL = 1e-14  # pressure computed from height above 86 km
TIME_ATOL = 1.2e-10  # days: one unit in the last place of a datenum


def _folder():
    try:
        return reference_dir()
    except (DataDirNotSetError, FileNotFoundError) as error:
        pytest.skip(f"no reference folder: {error}")


@pytest.fixture(scope="module")
def reference():
    file = _folder() / "get_limbsounders_reference.mat"
    if not file.is_file():
        pytest.skip(f"{file} not found")
    digest = hashlib.sha256(file.read_bytes()).hexdigest()
    assert digest == REFERENCE_SHA256, f"{file} is not the pinned reference"
    return read_mat(file)["Out"]


@pytest.fixture(scope="module")
def syn_dir():
    folder = _folder() / "get_limbsounders_synthetic"
    if not folder.is_dir():
        pytest.skip(f"{folder} not found")
    return folder


D1, D2 = "2007-01-01", "2007-01-02"


def synthetic_cases(syn):
    """The MATLAB test's synthetic calls, as (name, args, options)."""
    ml = {"misc_path": syn, "misc_identifier": "synlin"}
    mb = {"misc_path": syn, "misc_identifier": "synbad"}
    mp = {"hindley_pws": True, "misc_path": syn, "misc_identifier": "synpw"}
    mx = {"misc_path": syn, "misc_identifier": "synext"}
    return {
        "Syn_lin": ([D1, "Misc"], ml),
        "Syn_lin_2days_source": ([[D1, D2], "Misc"], {**ml, "file_source": True}),
        "Syn_bad_origz": ([D1, "Misc"], {**mb, "original_z": True}),
        "Syn_lin_origz_hr": ([D1, "Misc"],
                             {**ml, "original_z": True, "height_range": (20, 40)}),
        "Syn_pw_origz": ([D1, "Misc"], {**mp, "original_z": True}),
        "Syn_pw_origz_hr": ([D1, "Misc"],
                            {**mp, "original_z": True, "height_range": (20, 40)}),
        "Syn_pw": ([D2, "Misc"], mp),
        "Syn_3d": ([[D1, D2], "Misc"],
                   {"misc_path": syn, "misc_identifier": "syn3d",
                    "additional_vars": ["Extra3"]}),
        "Syn_ext": ([[D1, D2], "Misc"],
                    {**mx, "additional_vars": ["PerPointT", "PerProf", "PerLev"]}),
        "Syn_ext_ignore": ([D1, "Misc"], {**mx, "additional_vars": ["Meta", "TempMeta"],
                                          "ignore_non_modal": True}),
        "Syn_pw_ignore": ([D2, "Misc"], {**mp, "additional_vars": ["Meta"],
                                         "ignore_non_modal": True}),
        "Syn_flag": ([D1, "Misc"], {**mx, "additional_vars": ["FlagVar"]}),
        "Syn_lon_origz": ([D1, "Misc"], {"misc_path": syn, "misc_identifier": "synlon",
                                         "original_z": True}),
        "Syn_sampler": ([D1, "Sampler"], {"sample_path": Path(syn) / "sampler"}),
    }


SYNTHETIC = list(synthetic_cases(Path(".")))


@pytest.mark.parametrize("case", SYNTHETIC)
def test_synthetic_against_matlab(reference, syn_dir, case):
    args, options = synthetic_cases(syn_dir)[case]
    got, _ = _load_limb_arrays(*args, **_defaults(options))
    # Sampler files hold pressure, so only the others compute it from height
    tol = {} if case == "Syn_sampler" else {"Pres": ("rtol", PRESSURE_RTOL)}
    assert_same_as_matlab(got, reference[case], _interp_pressure_tol(tol, options))


# real data: (dataset, kind, args, options), as in the MATLAB test

HIRDLS_TOL = {"Lat": ("atol", POSITION_ATOL), "Lon": ("atol", POSITION_ATOL),
              "dx": ("rtol", DX_RTOL)}
REAL = {
    "HIRDLS": ("HIRDLS", "raw", [D1, "HIRDLS"], {"file_source": True}, HIRDLS_TOL),
    "HIRDLS_2days": ("HIRDLS", "raw", [[D1, D2], "HIRDLS"], {"lat_range": (-20, 20)},
                     HIRDLS_TOL),
    "HIRDLS_origz_hr": ("HIRDLS", "raw", [D1, "HIRDLS"],
                        {"original_z": True, "height_range": (20, 60)}, HIRDLS_TOL),
    "HIRDLS_scanup": ("HIRDLS", "raw", [D1, "HIRDLS"],
                      {"additional_vars": ["ScanUpFlag"], "lat_range": (-20, 20)},
                      HIRDLS_TOL),
    "HIRDLS_latlon_dx": ("HIRDLS", "raw", [D1, "HIRDLS"],
                         {"lat_range": (-30, 30), "lon_range": (-90, 90), "dx": True},
                         HIRDLS_TOL),
    "HIRDLS_pws": ("HIRDLS", "pws", [[D1, D2], "HIRDLS"],
                   {"hindley_pws": True, "lat_range": (-20, 20)},
                   {"Pres": ("rtol", PRESSURE_RTOL)}),
    "HIRDLS_extras": ("HIRDLS", "raw", [D1, "HIRDLS"],
                      {"additional_vars": ["TemperaturePrecision",
                                           "OrbitAscendingFlag"],
                       "lat_range": (-20, 20)}, HIRDLS_TOL),
    "SABER": ("SABER", "raw", [D1, "SABER"], {"file_source": True}, {}),
    "SABER_origz": ("SABER", "raw", [D1, "SABER"], {"original_z": True}, {}),
    "MLS": ("MLS", "raw", [D1, "MLS"], {"file_source": True}, {}),
    "MLS_extras": ("MLS", "raw", [D1, "MLS"],
                   {"additional_vars": ["L2gpPrecision", "Quality"],
                    "lat_range": (-20, 20)}, {}),
    "ACE": ("ACE", "raw", [D1, "ACE"], {"file_source": True}, {}),
    "GNSS": ("GNSS", "raw", [D1, "GNSS"], {"file_source": True}, {}),
    "GNSS_origz": ("GNSS", "raw", [D1, "GNSS"],
                   {"original_z": True, "lat_range": (-30, 30)}, {}),
    "MIPAS": ("MIPAS", "raw", [D2, "MIPAS"], {"file_source": True}, {}),
    "MIPAS_extras": ("MIPAS", "raw", [D2, "MIPAS"],
                     {"additional_vars": ["TEM_SD", "sza"]}, {}),
    "SOFIE": ("SOFIE", "raw", ["2007-06-01", "SOFIE"], {"file_source": True}, {}),
    # ALOMAR files hold no pressure, so it is computed from height
    "ALOMAR": ("ALOMAR", "raw", ["2007-01-23", "ALOMAR"], {"file_source": True},
               {"Pres": ("rtol", PRESSURE_RTOL)}),
    "ALOMAR_unc": ("ALOMAR", "raw", ["2007-01-23", "ALOMAR"],
                   {"additional_vars": ["T1_lidar_unc"]},
                   {"Pres": ("rtol", PRESSURE_RTOL)}),
    # added after the 2026-11-25 sweep (tests/matlab/sweep)
    "MLS_origz": ("MLS", "raw", [D1, "MLS"],
                  {"original_z": True, "lat_range": (-20, 20)},
                  {"Alt": ("rtol", HEIGHT_RTOL)}),
    "ALOMAR_origz_high": ("ALOMAR", "raw", ["2005-01-15", "ALOMAR"],
                          {"original_z": True},
                          {"Pres": ("rtol", HIGH_PRESSURE_RTOL),
                           "Time": ("atol", TIME_ATOL)}),
    "SABER_dx": ("SABER", "raw", [D1, "SABER"], {"dx": True, "lat_range": (-20, 20)},
                 {"dx": ("rtol", DX_RTOL)}),
    "MLS_2026_dup": ("MLS", "raw", ["2026-01-01", "MLS"],
                     {"file_source": True, "lat_range": (-20, 20)}, {}),
}


def test_every_reference_case_is_tested(reference):
    # AIRS and IASI are not ported; the MATLAB reference has no cases for them
    assert sorted(reference) == sorted(SYNTHETIC + list(REAL))


@pytest.mark.parametrize("case", list(REAL))
# (MLS_2026_dup: the warning about two files for the day is expected, and is
# tested in test_limb_known.py)
@pytest.mark.filterwarnings("ignore:.*more than one matching MLS file")
def test_real_data_against_matlab(reference, case):
    dataset, kind, args, options, tol = REAL[case]
    if not dataset_dir(dataset, kind, must_exist=False).is_dir():
        pytest.skip(f"no {dataset} {kind} data")
    want = reference[case]
    tol = _interp_pressure_tol(tol, options)
    if dataset in ("ACE", "GNSS"):
        got, want, tol = _python_changes(dataset, args, options, want, tol)
    else:
        if dataset in ("MLS", "SOFIE"):  # (Q15, Q17: MATLAB applies no screening)
            options = {**options, "instrument_qc": False}
        got, _ = _load_limb_arrays(*args, **_defaults(options))
        want = _q12_rows(want, options)
    assert_same_as_matlab(got, want, tol)


def _q12_rows(want, options):
    """The MATLAB rows a box keeps under Q12: at least one level whose
    latitude and longitude are both finite and inside the box (MATLAB kept
    profiles whose extremes overlapped it, possibly at different levels)."""
    if "lat_range" not in options and "lon_range" not in options:
        return want
    la0, la1 = sorted(options.get("lat_range", (-90, 90)))
    lo0, lo1 = options.get("lon_range", (-180, 180))
    lat, lon = np.asarray(want["Lat"], float), np.asarray(want["Lon"], float)
    with np.errstate(invalid="ignore"):
        inside = ((lat >= la0) & (lat <= la1) & (lon >= lo0) & (lon <= lo1))
    keep = np.flatnonzero(inside.any(axis=1))
    out = {k: (np.asarray(v)[keep] if k != "OriginalFiles" else v)
           for k, v in want.items()}
    if "dx" in out:  # (distance from the previous profile kept, where that changed)
        moved = np.flatnonzero(np.diff(keep) > 1) + 1
        p1, l1 = np.deg2rad(out["Lat"][moved - 1]), np.deg2rad(out["Lon"][moved - 1])
        p2, l2 = np.deg2rad(out["Lat"][moved]), np.deg2rad(out["Lon"][moved])
        h = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(
            (l2 - l1) / 2) ** 2
        h = np.clip(h, 0, 1)
        out["dx"] = np.array(out["dx"], dtype=float)
        out["dx"][moved] = 6371.0 * 2 * np.arctan2(np.sqrt(h), np.sqrt(1 - h))
    return out


def _python_changes(dataset, args, options, want, tol):
    """(got, want, tol) for ACE and GNSS: the MATLAB output with the
    deliberate differences applied (see the module docstring)."""
    want = {k: (np.array(v, dtype=float) if k not in ("OriginalFiles",) else v)
            for k, v in want.items()}
    if dataset == "ACE":
        got, _ = _load_limb_arrays(*args, **_defaults(options))
        want["Pres"] = want["Pres"] * 1013.25
        return got, want, {**tol, "Pres": tol.get("Pres", ("rtol", PRESSURE_RTOL))}

    # GNSS: which file row each profile is (Python: from its SourceProf)
    got, _ = _load_limb_arrays(*args, **_defaults({**options, "file_source": True}))
    got_rows = np.asarray(got["SourceProf"])[:, 0].astype(int) - 1
    if not options.get("file_source"):
        for name in ("SourceProf", "SourceFile", "OriginalFiles"):
            got.pop(name)
    day = args[0]
    year, doy = day[:4], (np.datetime64(day) - np.datetime64(day[:4] + "-01-01")
                          ).astype(int) + 1
    content = read_mat(dataset_dir("GNSS") / year / f"merged_ro_{year}_{doy:03d}.mat")
    meta = content["MetaData"]
    lat_file = np.asarray(content["Lat"], dtype=float)
    lon_file = np.asarray(content["Lon"], dtype=float)
    no_position = ~(np.isfinite(lat_file) & np.isfinite(lon_file)).any(axis=1)
    if "SourceProf" in want:
        want_rows = want["SourceProf"][:, 0].astype(int) - 1
    else:  # MATLAB's rows: the box test on the file's positions (NaN kept)
        want_rows = np.flatnonzero(_matlab_box(lat_file, lon_file, options))

    # rows left out by Python: repeated occultations, or (with a box) profiles
    # whose MetaData position is outside it
    occid = meta["occid"].astype(str).to_numpy()
    kept = set(got_rows.tolist())
    assert kept <= set(want_rows.tolist())
    for row in sorted(set(want_rows.tolist()) - kept):
        repeated = any(r in kept for r in np.flatnonzero(occid == occid[row]))
        outside = no_position[row] and not _matlab_box(
            np.array([[meta["latitude"][row]]], dtype=float),
            np.array([[meta["longitude"][row]]], dtype=float), options)[0]
        assert repeated or outside, f"file row {row} left out without reason"
    index = np.searchsorted(want_rows, got_rows)
    want = {k: (v[index] if k != "OriginalFiles" else v) for k, v in want.items()}

    want["Pres"] = want["Pres"] / 100
    # times: the occultation id's where the file's is hour-rounded
    file_time = np.asarray(meta["time"], dtype=float)
    stamps = [o[-12:] for o in occid]
    id_time = np.array([datenum(int(s[:4]), int(s[4:6]), int(s[6:8]), int(s[8:10]),
                                int(s[10:12])) for s in stamps])
    minutes = np.abs(file_time - id_time) * 1440
    on_hour = np.abs(file_time * 24 - np.round(file_time * 24)) < 1e-7
    off = (minutes > 10) | (on_hour & (minutes > 1))
    present = np.isfinite(want["Time"])  # (points the interpolation reaches)
    for i, row in enumerate(got_rows):
        if off[row]:
            want["Time"][i] = np.where(present[i], id_time[row], np.nan)
        if no_position[row]:
            want["Lat"][i] = np.where(present[i], meta["latitude"][row], np.nan)
            want["Lon"][i] = np.where(present[i], meta["longitude"][row], np.nan)
    return got, want, {**tol, "Lat": ("atol", POSITION_ATOL),
                       "Lon": ("atol", POSITION_ATOL),
                       "Pres": tol.get("Pres", ("rtol", PRESSURE_RTOL)),
                       "Time": ("atol", TIME_ATOL)}


def _matlab_box(lat, lon, options):
    """MATLAB's profile box test: a profile is kept unless its latitude or
    longitude extremes lie wholly outside the ranges (NaN-only profiles kept)."""
    la0, la1 = options.get("lat_range", (-90, 90))
    lo0, lo1 = options.get("lon_range", (-180, 180))
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        outside = ((np.nanmin(lon, axis=1) > max(lo0, lo1))
                   | (np.nanmax(lon, axis=1) < min(lo0, lo1))
                   | (np.nanmin(lat, axis=1) > max(la0, la1))
                   | (np.nanmax(lat, axis=1) < min(la0, la1)))
    return ~outside


def _interp_pressure_tol(tol, options):
    """tol, with the log-pressure tolerance for Pres when heights are interpolated."""
    if options.get("original_z"):
        return tol
    return {**tol, "Pres": ("rtol", LOG_PRESSURE_RTOL)}


def assert_same_as_matlab(got, want, tol=None):
    """Field-by-field comparison with a MATLAB output struct (as a dict).
    tol: {field: ("rtol" or "atol", value)}; other fields must be bit-identical."""
    tol = tol or {}
    assert sorted(got) == sorted(want), "different fields"
    for name in want:
        if name == "OriginalFiles":
            want_files = [str(np.asarray(f).squeeze())
                          for f in np.asarray(want[name], dtype=object).ravel()]
            assert list(got[name]) == want_files, name
            continue
        g = np.asarray(got[name])
        w = np.asarray(want[name])
        assert g.shape == w.shape, f"{name}: shape {g.shape} vs MATLAB {w.shape}"
        np.testing.assert_array_equal(np.isnan(g), np.isnan(w), err_msg=name)
        if name in tol:
            kind, value = tol[name]
            rtol, atol = (value, 0) if kind == "rtol" else (0, value)
            np.testing.assert_allclose(g, w, rtol=rtol, atol=atol, err_msg=name)
        else:
            np.testing.assert_array_equal(g, w, err_msg=name)


def _defaults(options):
    """Fill in the public function's defaults (the private layer needs them all)."""
    import inspect

    from atmosgw.curate.obs.limb import load_limb_sounders

    params = inspect.signature(load_limb_sounders).parameters
    full = {k: p.default for k, p in params.items()
            if p.default is not inspect.Parameter.empty}
    full.update(options)
    return full
