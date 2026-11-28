"""analysis.detrend and analysis.smooth against the group's MATLAB.

References (in reference_dir; SHA-256 pinned below, a different file fails):
- detrend_sgolay_ref.mat, made by tests/matlab/make_detrend_reference.m with
  MATLAB R2024b Update 6, corwin365/MatlabFunctions-Common 52518fc and
  -Personal 5112703: func_filter_sgolay (HIRDLS and SABER 2007-01-01, five
  length/order cases, SABER also on its original descending, irregular
  heights), func_filter_hindley23 (with and without H23_OutRem) and smoothn
  (box and Gaussian, with and without a NaN).
- prep_airs_3d_reference.mat (see tests/test_nadir_matlab.py): ret_temp, Tp
  and BG, [cross x along x level], for every AIRS and IASI case with Tp
  (the 4th-order cross-track polynomial of airs_4dp_detrend;
  AIRS_presmooth with PreSmooth [3,3,1]).

Tolerances:
- array level, float64: KELVIN = 1e-10 K, rounding only. Largest seen
  (2026-11-28): Savitzky-Golay 4.6e-13 K, polynomial 6.9e-13 K, smoothing
  1.8e-13 K, Hindley23 identical. (The polynomial design matrix of pixel
  numbers to the 4th power has a condition number near 1e9, so MATLAB's
  backslash and LAPACK's gelsy differ in the last bits.) NaN patterns must
  match exactly.
- Dataset level (outputs float32): profiles, whose inputs stay float64,
  PROFILE_F32 = 3.1e-5 K, one float32 spacing at 256-512 K (the rounding of
  the output is at most half of that); swaths, whose inputs are float32
  (the format), SWATH_F32 = 1e-4 K: the input rounding (up to 1.5e-5 K)
  carried through the fit plus the output rounding (largest seen 2.9e-5 K).

Differences from MATLAB tested here as such: Hindley23 background NaN where
every mode is NaN (MATLAB: 0; no such point in this reference).

Skipped (not passed) if reference_dir is not set or a file is absent.
"""

import hashlib
import importlib

import numpy as np
import pytest
from test_data_formats import example_swath
from test_detrend import example_profile, example_pw_profile

from atmosgw.curate.analysis import detrend, smooth_array
from atmosgw.curate.data import states, validate
from atmosgw.curate.paths import DataDirNotSetError, reference_dir
from atmosgw.curate.readers import read_mat

DT = importlib.import_module("atmosgw.curate.analysis.detrend")
SGOLAY_SHA256 = "03d0e9b639ad946a8606d06ed7023055e0453deda6618e9a8c4069b1d21c03ac"
AIRS_SHA256 = "a56e4f892f0791c47b789a1fdc4ccf5287401396505d1b772f83b94a4001725a"
KELVIN = 1e-10  # K, array level (see the module docstring)
PROFILE_F32 = 3.1e-5  # K, profile Datasets (float32 outputs)
SWATH_F32 = 1e-4  # K, swath Datasets (float32 inputs and outputs)
SG_CASES = ["SG_HIRDLS_18_2", "SG_SABER_18_2", "SG_SABER_origz_18_2",
            "SG_HIRDLS_5_2", "SG_HIRDLS_10_3"]
POLY_CASES = ["AIRS_default", "AIRS_dateline", "AIRS_midnight", "AIRS_dXT_dAT",
              "AIRS_nearest", "AIRS_daynight", "AIRS_presmooth",
              "AIRS_novertinterp", "AIRS_detrend1", "IASI_full", "IASI_midnight"]
CXL = ("cross_track", "along_track", "level")  # MATLAB's order
POLY = {"order": 4, "dim": "cross_track", "abscissa": "index", "min_points": None}


def _reference(name, sha256):
    try:
        file = reference_dir() / name
    except (DataDirNotSetError, FileNotFoundError) as error:
        pytest.skip(f"no reference folder: {error}")
    if not file.is_file():
        pytest.skip(f"{file} not found")
    digest = hashlib.sha256(file.read_bytes()).hexdigest()
    assert digest == sha256, f"{file} is not the pinned reference"
    return read_mat(file)["Ref"]


@pytest.fixture(scope="module")
def sgolay():
    return _reference("detrend_sgolay_ref.mat", SGOLAY_SHA256)


@pytest.fixture(scope="module")
def airs():
    out = _reference("prep_airs_3d_reference.mat", AIRS_SHA256)["Out"]
    return {name: case.get("Data", case) for name, case in out.items()}


def scalar(x):
    return np.ravel(x)[0]  # (read_mat returns MATLAB scalars as arrays)


def assert_close(got, want, atol, what):
    got, want = np.asarray(got, dtype=float), np.asarray(want, dtype=float)
    assert got.shape == want.shape, f"{what}: shape {got.shape} vs {want.shape}"
    np.testing.assert_array_equal(np.isnan(got), np.isnan(want), err_msg=what)
    np.testing.assert_allclose(got, want, rtol=0, atol=atol, equal_nan=True,
                               err_msg=what)


def sg_options(case):
    return {"order": int(scalar(case["SGOrder"])),
            "length_km": float(scalar(case["SGLength"])), "min_points": None}


def swath_dataset(ret_temp, ret_z):
    """A valid 3-D swath Dataset from MATLAB [cross x along x level]."""
    temp = np.transpose(np.asarray(ret_temp), (1, 0, 2)).astype(np.float32)
    ds = example_swath(*temp.shape)
    dims = ds["temperature"].dims
    flag = states.new_flags(temp.shape)
    states.set_state(flag, np.isnan(temp), "missing")
    ds["temperature"] = (dims, temp, ds["temperature"].attrs)
    ds["temperature_flag"] = (dims, flag, states.flag_attrs())
    ds["altitude"] = (("level",), np.ravel(ret_z).astype(np.float64),
                      {"units": "km"})
    return ds


def along_first(x):
    return np.transpose(np.asarray(x), (1, 0, 2))


# ------------------------------------------------------- array level


@pytest.mark.parametrize("name", SG_CASES)
def test_savgol_arrays(sgolay, name):
    case = sgolay[name]
    temp = np.asarray(case["Temp"], dtype=float)
    bg, _ = DT._savgol_step(case["Alt"], temp, ("profile", "level"),
                            sg_options(case), False)
    assert_close(bg, case["BG"], KELVIN, f"{name} BG")
    assert_close(temp - bg, case["Tp"], KELVIN, f"{name} Tp")


@pytest.mark.parametrize("name", ["H23_OutRem0", "H23_OutRem1"])
def test_hindley23_arrays(sgolay, name):
    given, want = sgolay["H23_In"], sgolay[name]
    waves = np.asarray(given["Temp_PW"], dtype=float)
    some = np.any(np.isfinite(waves), axis=-1)
    bg = np.where(some, np.nansum(waves, axis=-1), np.nan)  # (as detrend)
    # perturbation: the stored residual, exactly (the outlier removal of
    # H23_OutRem changes nothing here: the loader applied the same limits)
    np.testing.assert_array_equal(given["Temp_Residual"], want["Tp"])
    np.testing.assert_array_equal(bg[some], want["BG"][some])
    # the documented difference: MATLAB's nansum gives 0 where every mode is NaN
    np.testing.assert_array_equal(want["BG"][~some], 0)


SMOOTH_CASES = {
    "Box_3_5_1": ("X", (3, 5, 1), "box", 0.65),
    "Box_3_3_3": ("X", (3, 3, 3), "box", 0.65),
    "BoxNaN": ("XNaN", (3, 5, 1), "box", 0.65),
    "Gauss_5_3_3": ("X", (5, 3, 3), "gaussian", 0.65),
    "Gauss_5_5_1_std": ("X", (5, 5, 1), "gaussian", "Gauss_std"),
}


@pytest.mark.parametrize("name", list(SMOOTH_CASES))
def test_smoothing_arrays(sgolay, name):
    ref = sgolay["Smooth"]
    source, size, kernel, sigma = SMOOTH_CASES[name]
    if isinstance(sigma, str):
        sigma = list(np.ravel(ref[sigma]))
    got = smooth_array(ref[source], size, kernel, sigma)
    assert_close(got, ref[name], KELVIN, name)


@pytest.mark.parametrize("name", POLY_CASES)
def test_polynomial_arrays(airs, name):
    case = airs[name]
    assert "Tp" in case and "BG" in case, f"{name} has no Tp/BG"
    temp = np.asarray(case["ret_temp"], dtype=float)
    bg, _ = DT._polynomial_step(None, temp, CXL, POLY, False)
    assert_close(bg, case["BG"], KELVIN, f"{name} BG")
    assert_close(temp - bg, case["Tp"], KELVIN, f"{name} Tp")


def test_every_case_with_tp_is_compared(airs):
    assert sorted(n for n, c in airs.items() if "Tp" in c) == sorted(POLY_CASES)


def test_presmooth_is_smoothing_then_polynomial(airs):
    # AIRS_presmooth's ret_temp is the smoothed field (PreSmooth [3,3,1] in
    # test_prep_airs_3d.m): smoothn of AIRS_default's ret_temp
    default, smoothed = airs["AIRS_default"], airs["AIRS_presmooth"]
    assert not np.array_equal(default["ret_temp"], smoothed["ret_temp"])
    got = smooth_array(default["ret_temp"], (3, 3, 1))
    assert_close(got, smoothed["ret_temp"], KELVIN, "smoothed ret_temp")
    bg, _ = DT._polynomial_step(None, got, CXL, POLY, False)
    assert_close(bg, smoothed["BG"], KELVIN, "presmooth BG")
    assert_close(got - bg, smoothed["Tp"], KELVIN, "presmooth Tp")


# ----------------------------------------------------- Dataset level


@pytest.mark.parametrize("name", SG_CASES)
def test_savgol_dataset(sgolay, name):
    case = sgolay[name]
    ds = example_profile(case["Temp"], case["Alt"])
    options = sg_options(case)
    out = detrend(ds, method="savgol", order=options["order"],
                  length_km=options["length_km"], mask_tropopause_kink=False)
    assert_close(out["temperature_background"], case["BG"], PROFILE_F32, "BG")
    assert_close(out["temperature_perturbation"], case["Tp"], PROFILE_F32, "Tp")
    validate(out, "profile")


def test_hindley23_dataset(sgolay):
    given, want = sgolay["H23_In"], sgolay["H23_OutRem0"]
    ds = example_pw_profile(given["Temp_Residual"], given["Temp_PW"], given["Alt"])
    out = detrend(ds)  # (the default for planetary-wave data)
    assert_close(out["temperature_perturbation"], want["Tp"], PROFILE_F32, "Tp")
    assert_close(out["temperature_background"], want["BG"], PROFILE_F32, "BG")
    validate(out, "profile")


def test_polynomial_dataset(airs):
    case = airs["AIRS_default"]
    out = detrend(swath_dataset(case["ret_temp"], case["ret_z"]))
    assert_close(out["temperature_background"], along_first(case["BG"]),
                 SWATH_F32, "BG")
    assert_close(out["temperature_perturbation"], along_first(case["Tp"]),
                 SWATH_F32, "Tp")
    validate(out, "swath")


def test_smooth_then_polynomial_dataset(airs):
    default, smoothed = airs["AIRS_default"], airs["AIRS_presmooth"]
    out = detrend(swath_dataset(default["ret_temp"], default["ret_z"]),
                  steps=[("smooth", {"size": {"cross_track": 3,
                                              "along_track": 3}}),
                         "polynomial"])
    assert_close(out["temperature_background"], along_first(smoothed["BG"]),
                 SWATH_F32, "BG")
    assert_close(out["temperature_perturbation"], along_first(smoothed["Tp"]),
                 SWATH_F32, "Tp")
    validate(out, "swath")
