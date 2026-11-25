"""Tests for atmosgw.curate.solar.

MATLAB outputs: tests/data/matlab_swath_fixtures.mat (see
tests/matlab/make_swath_fixtures.m): ``which_airs_retrieval`` (the day/night
mask) and the solar zenith angle by the same formulas, on a 10 x 20 degree
grid at six times.

Tolerances against MATLAB: the angle within 1e-6 degrees (MATLAB works from
datenums of about 7.3e5, whose rounding, about 1e-10 days, moves the hour
angle by up to about 4e-8 degrees; the datetime64 route here is exact). The
mask must match exactly, except within 1e-6 degrees of the 96-degree
threshold, where that rounding can decide it.
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from atmosgw.curate.matlab_time import datenum_to_datetime
from atmosgw.curate.solar import day_retrieval, solar_zenith_angle

FIXTURES = Path(__file__).parent / "data" / "matlab_swath_fixtures.mat"


@pytest.fixture(scope="module")
def matlab():
    return loadmat(FIXTURES, squeeze_me=True, struct_as_record=False)


@pytest.mark.parametrize("use_datenum", [False, True])
def test_against_matlab(matlab, use_datenum):
    lat, lon = matlab["SunLat"], matlab["SunLon"]
    for k, t in enumerate(np.atleast_1d(matlab["SunTimes"])):
        time = t if use_datenum else datenum_to_datetime(t, "us")
        sza = solar_zenith_angle(time, lat, lon, datenum=use_datenum)
        want = matlab["SunSZA"][:, :, k]
        np.testing.assert_allclose(sza, want, rtol=0, atol=1e-6)
        mask = day_retrieval(time, lat, lon, datenum=use_datenum)
        decided = np.abs(want - 96) > 1e-6
        np.testing.assert_array_equal(mask[decided],
                                      matlab["SunMask"][:, :, k][decided])


def test_sun_overhead_at_june_solstice_noon():
    # 2020-06-21 12:00 UTC: the Sun is near the Tropic of Cancer at about 0 E
    # (equation of time about -1.6 min, i.e. 0.4 degrees of longitude)
    t = np.datetime64("2020-06-21T12:00")
    assert solar_zenith_angle(t, 23.44, 0.0) < 0.5
    assert abs(solar_zenith_angle(t, -66.56, 0.0) - 90) < 0.5  # on the terminator


def test_antipode_is_supplement():
    t = np.datetime64("2011-03-01T05:17")
    lat, lon = 37.0, -122.0
    a = solar_zenith_angle(t, lat, lon)
    b = solar_zenith_angle(t, -lat, lon + 180)
    assert abs(a + b - 180) < 1e-9


def test_mask_values_and_shape():
    t = np.datetime64("2007-06-21T00:00")
    lat, lon = np.meshgrid(np.arange(-80, 81, 40.0), np.arange(-180, 180, 90.0))
    mask = day_retrieval(t, lat, lon)
    assert mask.shape == lat.shape
    assert set(np.unique(mask)) <= {0.0, 1.0}
    assert np.isnan(day_retrieval(t, np.nan, 0.0))
