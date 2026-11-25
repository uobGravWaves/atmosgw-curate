"""Tests for atmosgw.curate.atmosphere (U.S. Standard Atmosphere 1976, [R1]).

MATLAB outputs: tests/data/matlab_compat_fixtures.mat (see
tests/matlab/make_compat_fixtures.m), from the rewritten h2p.m and p2h.m
(corwin365/MatlabFunctions, see the fixture's Meta), MATLAB R2024b Update 6.

Tolerances against MATLAB: pressure 1e-14 relative, height 1e-12 km absolute.
The formulas and order of operations are the same, but power, exp and log
round differently in the last bits between MATLAB and NumPy, and the table
interpolation is done in the opposite direction: on 2026-11-25, over 4000
values, the largest differences were 2.0e-15 relative (pressure) and 5.7e-14
km (height). The NaN pattern must match exactly.

Known answers: pressures printed in Table I of [R1] (5 digits), so 1e-4
relative.
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from atmosgw.curate.atmosphere import (
    BASE_HEIGHTS,
    BASE_PRESSURES,
    EARTH_RADIUS,
    TABLE_HEIGHTS,
    TABLE_PRESSURES,
    height_to_pressure,
    pressure_to_height,
)

FIXTURES = Path(__file__).parent / "data" / "matlab_compat_fixtures.mat"
PRESSURE_RTOL = 1e-14
HEIGHT_ATOL = 1e-12  # km

# [R1] Table I (geometric altitude), pp. 51, 65, 67, 68: Z (km) -> P (hPa)
TABLE_I = {-5: 1.7776e3, 0: 1013.25, 40: 2.8714, 50: 7.9779e-1, 55: 4.2525e-1,
           60: 2.1958e-1, 70: 5.2209e-2, 80: 1.0524e-2, 85: 4.4568e-3,
           85.5: 4.0802e-3, 86: 3.7338e-3, 100: 3.2011e-4, 149: 4.7421e-6}


@pytest.fixture(scope="module")
def matlab():
    return loadmat(FIXTURES, squeeze_me=True, struct_as_record=False)


def assert_matches(got, want, rtol=0, atol=0):
    got, want = np.asarray(got), np.asarray(want, dtype=float)
    assert got.shape == want.shape
    np.testing.assert_array_equal(np.isnan(got), np.isnan(want))
    np.testing.assert_allclose(got, want, rtol=rtol, atol=atol)


def test_height_to_pressure_against_matlab(matlab):
    assert_matches(height_to_pressure(matlab["H2pIn"]), matlab["H2pOut"],
                   rtol=PRESSURE_RTOL)
    assert_matches(height_to_pressure(matlab["H2pMatrix"]), matlab["H2pMatrixOut"],
                   rtol=PRESSURE_RTOL)


def test_pressure_to_height_against_matlab(matlab):
    assert_matches(pressure_to_height(matlab["P2hIn"]), matlab["P2hOut"],
                   atol=HEIGHT_ATOL)
    assert_matches(pressure_to_height(matlab["P2hMatrix"]), matlab["P2hMatrixOut"],
                   atol=HEIGHT_ATOL)


# known answers from [R1]


@pytest.mark.parametrize("z", list(TABLE_I))
def test_table_i_pressures(z):
    assert height_to_pressure(z) == pytest.approx(TABLE_I[z], rel=1e-4)


def test_layer_base_pressures():
    # base pressures computed from eqs (33a)/(33b) against Table I's values at
    # the geopotential base heights (geopotential pages, 5 digits)
    printed = [1013.25, 226.32, 54.749, 8.6802, 1.1091, 0.66939, 0.039564]
    np.testing.assert_allclose(BASE_PRESSURES, printed, rtol=1e-4)
    z = EARTH_RADIUS * BASE_HEIGHTS[:-1] / (EARTH_RADIUS - BASE_HEIGHTS[:-1])
    np.testing.assert_allclose(height_to_pressure(z), BASE_PRESSURES, rtol=1e-13)


def test_table_values_reproduced_exactly():
    np.testing.assert_allclose(height_to_pressure(TABLE_HEIGHTS), TABLE_PRESSURES,
                               rtol=1e-15)


def test_continuous_at_86_km():
    below, above = height_to_pressure([86 - 1e-6, 86.0])
    assert below == pytest.approx(above, rel=2e-6)  # equations vs table: 1.2e-6


def test_range_and_invalid_input():
    assert np.all(np.isnan(height_to_pressure([-5.001, 149.001, np.nan, np.inf])))
    assert np.all(np.isfinite(height_to_pressure([-5.0, 149.0])))
    assert np.all(np.isnan(pressure_to_height([1800.0, 4.7e-6, 0.0, -1.0, np.nan])))
    assert pressure_to_height(1013.25) == 0.0
    assert pressure_to_height(1050.0) < 0  # below sea level is supported


def test_round_trip():
    z = np.linspace(-5, 149, 20001)
    np.testing.assert_allclose(pressure_to_height(height_to_pressure(z)), z,
                               rtol=0, atol=1e-12)


def test_monotonic():
    p = height_to_pressure(np.linspace(-5, 149, 20001))
    assert np.all(np.diff(p) < 0)


def test_scalar_shape_and_float32():
    assert isinstance(height_to_pressure(5.0), float)
    assert np.isnan(pressure_to_height(np.nan))
    assert height_to_pressure(np.ones((2, 3))).shape == (2, 3)
    # float32 input is computed in double, like any other input
    got = pressure_to_height(np.float32(0.1))
    assert got.dtype == np.float64
    assert got == pressure_to_height(float(np.float32(0.1)))
