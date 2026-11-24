"""Tests for atmosgw.curate.atmosphere (ports of h2p.m and p2h.m).

MATLAB outputs: tests/data/matlab_compat_fixtures.mat (see
tests/matlab/make_compat_fixtures.m), from h2p.m and p2h.m at
corwin365/MatlabFunctions commit c943399, MATLAB R2024b Update 6.

Tolerance against MATLAB: 1e-15 relative (a few units in the last place). The
formulas and order of operations are the same as MATLAB's, but exp() and log()
round differently in the last bit between MATLAB and NumPy: on 2026-11-24, 17
of 277 values differed, by at most 2.3e-16 relative. The NaN pattern must
match exactly.
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from atmosgw.curate.atmosphere import (
    LEVEL_HEIGHTS,
    LEVEL_PRESSURES,
    height_to_pressure,
    pressure_to_height,
)

FIXTURES = Path(__file__).parent / "data" / "matlab_compat_fixtures.mat"
RTOL = 1e-15


@pytest.fixture(scope="module")
def matlab():
    return loadmat(FIXTURES, squeeze_me=True, struct_as_record=False)


def assert_matches(got, want):
    got, want = np.asarray(got), np.asarray(want, dtype=float)
    assert got.shape == want.shape
    np.testing.assert_array_equal(np.isnan(got), np.isnan(want))
    np.testing.assert_allclose(got, want, rtol=RTOL, atol=0)


def test_height_to_pressure_against_matlab(matlab):
    assert_matches(height_to_pressure(matlab["H2pIn"]), matlab["H2pOut"])
    assert_matches(height_to_pressure(matlab["H2pMatrix"]), matlab["H2pMatrixOut"])


def test_pressure_to_height_against_matlab(matlab):
    # includes Q3: 1050 hPa after 500 hPa reuses the layer of the 500 hPa element
    assert_matches(pressure_to_height(matlab["P2hIn"]), matlab["P2hOut"])
    assert_matches(pressure_to_height(matlab["P2hMatrix"]), matlab["P2hMatrixOut"])


def test_pressure_to_height_single_precision_against_matlab():
    # float32 input: MATLAB p2h computes in single precision (fixture:
    # tests/matlab/make_p2h_single_fixture.m, MATLAB R2024b Update 6).
    # Tolerance: 1 float32 unit in the last place, for at most 3% of values.
    # Every step matches MATLAB except its single-precision log, which is not
    # correctly rounded and so cannot be reproduced exactly (on 2026-11-24, 7
    # of 355 heights differed by 1 unit; the 55 MLS levels all match, which
    # tests/test_limb_matlab.py checks bit for bit).
    f = loadmat(Path(__file__).parent / "data" / "matlab_p2h_single.mat",
                squeeze_me=True)
    assert f["P"].dtype == np.float32
    got = pressure_to_height(f["P"])
    assert got.dtype == np.float64
    ulp = np.spacing(np.abs(f["Z"]).astype(np.float32)).astype(np.float64)
    assert np.all(np.abs(got - f["Z"]) <= ulp)
    assert np.mean(got != f["Z"]) <= 0.03


def test_q3_first_element_above_1013_errors_like_matlab(matlab):
    assert matlab["P2hFirstHighErrors"]
    with pytest.raises(ValueError, match="1013.25"):
        pressure_to_height(np.array([1050.0, 500.0]))


def test_q3_previous_element_in_column_major_order():
    # column-major order is 500, 1050 (first column), then 5, 1 ...
    p = np.array([[500.0, 5.0], [1050.0, 1.0]])
    got = pressure_to_height(p)
    assert got[1, 0] == pressure_to_height(np.array([500.0, 1050.0]))[1]


def test_q2_negative_altitude_is_nan():
    assert np.isnan(height_to_pressure(-0.5))


# known answers


def test_level_pressures_exact():
    np.testing.assert_array_equal(height_to_pressure(LEVEL_HEIGHTS), LEVEL_PRESSURES)


def test_round_trip():
    z = np.linspace(0.1, 120, 500)
    np.testing.assert_allclose(pressure_to_height(height_to_pressure(z)), z, rtol=1e-12)


def test_scalar_and_nan():
    assert isinstance(height_to_pressure(5.0), float)
    assert np.isnan(pressure_to_height(np.nan))
