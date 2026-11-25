"""Regression tests: load_nadir_sounders (AIRS) against the MATLAB prep_airs_3d.

Reference: <reference_dir>/prep_airs_3d_reference.mat, made 2026-11-25 by
corwin365/MatlabFunctions-Tests test_prep_airs_3d.m with MATLAB R2024b
Update 6, from corwin365/MatlabFunctions@059037e (AIRS outputs unchanged at
f1a78c0, checked in MATLAB on 2026-11-25). Its SHA-256 is pinned below; a
different file fails the test rather than being used.

What is compared: the private layer (_load_nadir_arrays), whose output has
MATLAB-like names, [cross x along (x level)] arrays and datenum times, against
each AIRS case of the MATLAB test (2007-01-01: granule 53 with default,
spacing, nearest, no vertical interpolation, day/night and load-only options;
granule 3, which crosses the dateline; granules 240 and 241 across midnight).
MATLAB has no outlier screen for AIRS, so keep_outliers=True. Detrended
fields (Tp, BG) are compared when analysis.detrend is ported.

Tolerances (trig functions round differently in the last bits between MATLAB
and NumPy; see tests/test_swath.py):
- latitude and longitude: 1e-9 degrees (largest seen 2.6e-11);
- temperature: 1e-9 K (largest seen 8.8e-13);
- time: 1e-9 days, 0.1 ms (largest seen 3.5e-10);
- spacing: 1e-9 km;
- load-only, day/night, heights: identical.
NaN patterns and shapes must match exactly.

Skipped (not passed) if reference_dir is not set, the reference is absent, or
the AIRS data are not found.
"""

import hashlib
import inspect

import numpy as np
import pytest

from atmosgw.curate.obs.nadir import _load_nadir_arrays, load_nadir_sounders
from atmosgw.curate.paths import DataDirNotSetError, dataset_dir, reference_dir
from atmosgw.curate.readers import read_mat

REFERENCE_SHA256 = "a56e4f892f0791c47b789a1fdc4ccf5287401396505d1b772f83b94a4001725a"
DEG, KELVIN, DAYS, KM = 1e-9, 1e-9, 1e-9, 1e-9
D, D2 = "2007-01-01", "2007-01-02"
DEFAULTS = {k: p.default for k, p in inspect.signature(load_nadir_sounders)
            .parameters.items() if p.default is not inspect.Parameter.empty}
CASES = {
    "AIRS_default": dict(granules=[(D, 53)]),
    "AIRS_dateline": dict(granules=[(D, 3)]),
    "AIRS_midnight": dict(granules=[(D, 240), (D2, 1)]),
    "AIRS_dXT_dAT": dict(granules=[(D, 53)], dx_cross=10, dx_along=10),
    "AIRS_nearest": dict(granules=[(D, 53)], interp_type="nearest"),
    "AIRS_novertinterp": dict(granules=[(D, 53)], original_z=True),
    "AIRS_daynight": dict(granules=[(D, 53)], day_night=True),
    "AIRS_loadonly": dict(granules=[(D, 53)], original_grid=True),
}


@pytest.fixture(scope="module")
def reference():
    try:
        file = reference_dir() / "prep_airs_3d_reference.mat"
        airs = dataset_dir("AIRS")
    except (DataDirNotSetError, FileNotFoundError) as error:
        pytest.skip(f"no reference folder or AIRS data: {error}")
    if not file.is_file():
        pytest.skip(f"{file} not found")
    if not (airs / "2007" / "airs_2007d001.nc").is_file():
        pytest.skip("AIRS 2007-01-01 not found")
    digest = hashlib.sha256(file.read_bytes()).hexdigest()
    assert digest == REFERENCE_SHA256, f"{file} is not the pinned reference"
    return read_mat(file)["Ref"]["Out"]


def arrays(**options):
    data, _ = _load_nadir_arrays(None, "AIRS",
                                 **{**DEFAULTS, **options, "keep_outliers": True})
    return data


def assert_close(got, want, atol, what):
    got, want = np.asarray(got, dtype=float), np.asarray(want, dtype=float)
    assert got.shape == want.shape, f"{what}: shape {got.shape} vs {want.shape}"
    np.testing.assert_array_equal(np.isnan(got), np.isnan(want), err_msg=what)
    np.testing.assert_allclose(got, want, rtol=0, atol=atol, equal_nan=True,
                               err_msg=what)


@pytest.mark.parametrize("case", list(CASES))
def test_against_matlab(reference, case):
    ref = reference[case]
    want = ref.get("Data", ref)
    data = arrays(**CASES[case])
    exact = case == "AIRS_loadonly"
    temp = want["ret_temp"]
    if exact:  # MATLAB LoadOnly returns [level x cross x along]
        temp = np.transpose(temp, (1, 2, 0))
    assert_close(data["Lat"], want["l1_lat"], 0 if exact else DEG, "latitude")
    dlon = np.mod(data["Lon"] - want["l1_lon"] + 180, 360) - 180
    assert_close(dlon, np.zeros_like(dlon), 0 if exact else DEG, "longitude")
    assert_close(data["Time"], want["l1_time"], 0 if exact else DAYS, "time")
    assert_close(data["Temp"], temp, 0 if exact else KELVIN, "temperature")
    np.testing.assert_array_equal(np.ravel(data["Alt"]), np.ravel(want["ret_z"]))
    if "Spacing" in ref and "Data" in ref:
        assert_close(data["Spacing"], np.ravel(ref["Spacing"]), KM, "spacing")
    if case == "AIRS_daynight":
        np.testing.assert_array_equal(data["DayNight"], want["DayNightFlag"])
