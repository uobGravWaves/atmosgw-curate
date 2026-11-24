"""Tests for atmosgw.curate.matlab_compat, against MATLAB R2024b output.

The MATLAB outputs are in tests/data/matlab_compat_fixtures.mat, made by
tests/matlab/make_compat_fixtures.m (inputs are stored with the outputs).

Tolerances:
- linear, nearest, next, previous, unwrap, wrap_to_pi: bit-identical, including
  where the NaNs are.
- pchip, spline, makima: 1e-12 relative to the largest sample value. These use
  scipy, whose arithmetic differs from MATLAB's in the last bits; the NaN
  pattern must still match exactly.
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from atmosgw.curate.matlab_compat import interp1, unwrap, wrap_to_pi

FIXTURES = Path(__file__).parent / "data" / "matlab_compat_fixtures.mat"
SMOOTH_RTOL = 1e-12


@pytest.fixture(scope="module")
def matlab():
    return loadmat(FIXTURES, squeeze_me=True, struct_as_record=False)


def as_float(a):
    return np.atleast_1d(np.asarray(a, dtype=float))


def interp_cases():
    """(case name, method, extrapolate) for every stored MATLAB result."""
    f = loadmat(FIXTURES, squeeze_me=True, struct_as_record=False)
    return [
        (name, method, extrap)
        for name in f["Interp"]._fieldnames
        for method in f["Methods"]
        for extrap in (False, True)
    ]


@pytest.mark.parametrize("name,method,extrap", interp_cases())
def test_interp1_against_matlab(matlab, name, method, extrap):
    case = getattr(matlab["Interp"], name)
    want = as_float(getattr(case, method + ("_extrap" if extrap else "")))
    # MATLAB's default extrapolation is the case without 'extrap'
    got = interp1(case.X, case.V, case.XQ, method, extrapolate=True if extrap else None)
    got = np.atleast_1d(got)
    np.testing.assert_array_equal(np.isnan(got), np.isnan(want))
    if method in ("pchip", "spline", "makima"):
        scale = np.nanmax(np.abs(as_float(case.V)))
        np.testing.assert_allclose(got, want, rtol=0, atol=SMOOTH_RTOL * scale)
    else:
        np.testing.assert_array_equal(got, want)


@pytest.mark.parametrize("name", ["random", "exact", "nearpi", "dateline"])
def test_unwrap_against_matlab(matlab, name):
    got = unwrap(getattr(matlab["Unwrap"], name))
    np.testing.assert_array_equal(got, as_float(getattr(matlab["UnwrapOut"], name)))


def test_wrap_to_pi_against_matlab(matlab):
    np.testing.assert_array_equal(wrap_to_pi(matlab["WrapIn"]), matlab["WrapOut"])


def test_deg2rad_rad2deg_need_no_copies(matlab):
    np.testing.assert_array_equal(np.deg2rad(matlab["DegIn"]), matlab["DegToRad"])
    np.testing.assert_array_equal(np.rad2deg(matlab["RadIn"]), matlab["RadToDeg"])


# known answers (no MATLAB needed)


def test_interp1_linear_data_exact_and_nan_outside():
    x = np.array([5.0, 1.0, 3.0, 2.0])  # unsorted
    got = interp1(x, 2 * x + 1, np.array([0.0, 1.0, 1.5, 4.9, 5.0, 6.0]))
    np.testing.assert_allclose(got, [np.nan, 3.0, 4.0, 10.8, 11.0, np.nan], rtol=1e-15)


def test_interp1_descending_same_as_ascending():
    x = np.array([1.0, 2.0, 4.0, 7.0])
    v = np.array([3.0, -1.0, 2.0, 5.0])
    xq = np.linspace(0, 8, 33)
    methods = ("linear", "nearest", "next", "previous", "pchip", "spline", "makima")
    for method in methods:
        np.testing.assert_array_equal(
            interp1(x[::-1], v[::-1], xq, method), interp1(x, v, xq, method)
        )


def test_interp1_linear_between_equal_values_is_exact():
    # MATLAB returns exactly the value; (1-t)*v0 + t*v1 alone can be a bit off
    # (checked against MATLAB R2024b on constant and step-wise data, 2026-11-24)
    x = np.array([10.3, 12.3, 14.3, 16.3])
    v = np.full(4, -8.57142857142857)
    got = interp1(x, v, np.linspace(10.3, 16.3, 1001))
    assert np.all(got == v[0])


def test_interp1_nearest_tie_goes_up():
    x, v = np.array([0.0, 1.0]), np.array([10.0, 20.0])
    assert interp1(x, v, np.array([0.5]), "nearest")[0] == 20.0


def test_interp1_keeps_shape_and_nan_queries():
    xq = np.array([[1.5, np.nan], [2.5, 3.0]])
    got = interp1(np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0]), xq)
    assert got.shape == (2, 2)
    np.testing.assert_array_equal(got, [[1.5, np.nan], [2.5, 3.0]])


def test_interp1_rejects():
    x, v = np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0])
    with pytest.raises(NotImplementedError):
        interp1(x, v, x, "cubic")
    with pytest.raises(ValueError, match="distinct"):
        interp1(np.array([1.0, 2.0, 2.0]), v, x)
    with pytest.raises(ValueError, match="Unknown"):
        interp1(x, v, x, "quadratic")


def test_unwrap_removes_dateline_jump():
    lon = np.deg2rad(np.array([178.0, 179.0, -179.0, -178.0]))
    np.testing.assert_allclose(
        np.rad2deg(unwrap(lon)), [178.0, 179.0, 181.0, 182.0], rtol=1e-13
    )


def test_wrap_to_pi_odd_multiples():
    np.testing.assert_array_equal(
        wrap_to_pi(np.array([np.pi, -np.pi, 3 * np.pi, -3 * np.pi, 2 * np.pi])),
        [np.pi, -np.pi, np.pi, -np.pi, 0.0],
    )
    assert wrap_to_pi(1.0) == 1.0
