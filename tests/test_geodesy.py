"""Tests for atmosgw.curate.geodesy, against MATLAB R2024b output.

MATLAB outputs: tests/data/matlab_compat_fixtures.mat (see
tests/matlab/make_compat_fixtures.m): Mapping Toolbox azimuth, reckon and
km2deg with their spherical defaults, and nph_haversine.m
(corwin365/MatlabFunctions commit c943399).

Tolerances against MATLAB (trig functions round differently in the last bits
between MATLAB and NumPy; formulas are the standard ones, so small
differences in how MATLAB arranges them are also possible):
- km2deg: bit-identical.
- haversine: 1e-14 relative.
- azimuth, and reckon latitude and longitude: 1e-11 degrees absolute (about
  1 micrometre on the ground).
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from atmosgw.curate.geodesy import azimuth, haversine, km2deg, reckon

FIXTURES = Path(__file__).parent / "data" / "matlab_compat_fixtures.mat"
DEG_ATOL = 1e-11


@pytest.fixture(scope="module")
def matlab():
    return loadmat(FIXTURES, squeeze_me=True, struct_as_record=False)


def test_haversine_against_matlab(matlab):
    a, b = matlab["GeoA"], matlab["GeoB"]
    got = haversine(a[:, 0], a[:, 1], b[:, 0], b[:, 1])
    np.testing.assert_allclose(got, matlab["Haversine"], rtol=1e-14, atol=0)


def test_azimuth_against_matlab(matlab):
    a, b = matlab["GeoA"], matlab["GeoB"]
    got = azimuth(a[:, 0], a[:, 1], b[:, 0], b[:, 1])
    np.testing.assert_allclose(got, matlab["Azimuth"], rtol=0, atol=DEG_ATOL)
    # special cases (identical, antipodal, dateline, poles, cardinal) exactly
    np.testing.assert_array_equal(got[-11:][[0, 3, 4, 5, 6, 7, 8, 9, 10]],
                                  matlab["Azimuth"][-11:][[0, 3, 4, 5, 6, 7, 8, 9, 10]])


def test_reckon_against_matlab(matlab):
    m = matlab
    lat, lon = reckon(m["ReckonLat"], m["ReckonLon"], m["ReckonArc"], m["ReckonAz"])
    np.testing.assert_allclose(lat, matlab["ReckonLatOut"], rtol=0, atol=DEG_ATOL)
    np.testing.assert_allclose(lon, matlab["ReckonLonOut"], rtol=0, atol=DEG_ATOL)

    m = matlab["ReckonMatrix"]
    lat, lon = reckon(m.Lat, m.Lon, m.Arc, m.Az)
    assert lat.shape == (3, 4)
    np.testing.assert_allclose(lat, m.LatOut, rtol=0, atol=DEG_ATOL)
    np.testing.assert_allclose(lon, m.LonOut, rtol=0, atol=DEG_ATOL)


def test_km2deg_against_matlab(matlab):
    np.testing.assert_array_equal(km2deg(matlab["KmIn"]), matlab["KmToDeg"])


# known answers


def test_quarter_circumference():
    assert haversine(0, 0, 0, 90) == pytest.approx(np.pi * 6371 / 2, rel=1e-15)
    assert haversine(0, 0, 90, 0) == pytest.approx(np.pi * 6371 / 2, rel=1e-15)


def test_cardinal_directions():
    got = azimuth(0, 0, np.array([1, 0, -1, 0]), np.array([0, 1, 0, -1]))
    np.testing.assert_allclose(got, [0, 90, 180, 270], atol=1e-12)


def test_reckon_round_trip():
    rng = np.random.default_rng(1)
    lat = rng.uniform(-80, 80, 200)
    lon = rng.uniform(-180, 180, 200)
    arc = rng.uniform(0.01, 60, 200)
    az = rng.uniform(0, 360, 200)
    lat2, lon2 = reckon(lat, lon, arc, az)
    np.testing.assert_allclose(km2deg(haversine(lat, lon, lat2, lon2)), arc, atol=1e-9)
    back = np.mod(azimuth(lat, lon, lat2, lon2) - az + 180, 360) - 180
    np.testing.assert_allclose(back, 0, atol=1e-8)


def test_negative_arc_goes_backwards():
    lat_f, lon_f = reckon(10, 20, 1, 90)
    lat_b, lon_b = reckon(10, 20, -1, 270)
    assert lat_f == pytest.approx(lat_b, abs=1e-12)
    assert lon_f == pytest.approx(lon_b, abs=1e-12)


def test_haversine_matches_vector_method_and_matlab_example():
    # independent method: angle between unit vectors, atan2(|a x b|, a . b)
    rng = np.random.default_rng(2)
    lat1, lat2 = np.rad2deg(np.arcsin(rng.uniform(-1, 1, (2, 10000))))
    lon1, lon2 = rng.uniform(-180, 180, (2, 10000))
    p1, l1, p2, l2 = (np.deg2rad(a) for a in (lat1, lon1, lat2, lon2))
    a = np.stack([np.cos(p1) * np.cos(l1), np.cos(p1) * np.sin(l1), np.sin(p1)])
    b = np.stack([np.cos(p2) * np.cos(l2), np.cos(p2) * np.sin(l2), np.sin(p2)])
    angle = np.arctan2(np.linalg.norm(np.cross(a, b, axis=0), axis=0), (a * b).sum(0))
    np.testing.assert_allclose(
        haversine(lat1, lon1, lat2, lon2), 6371 * angle, rtol=0, atol=1e-9
    )
    # worked example in the nph_haversine.m header: ~11436.398 km
    assert haversine(63.2, -30, 13.95, 150) == pytest.approx(11436.398, abs=5e-4)


def test_q4_near_antipodal_gives_half_circle():
    # Q4 (fixed): rounding pushed a above 1 for some nearly antipodal pairs,
    # giving NaN; a is now clipped, so they give half a great circle.
    # Tolerance 1e-3 km: the formula is ill-conditioned near antipodes (a
    # rounding error of ~1e-16 in a becomes ~1e-8 rad, about 0.1 m), plus the
    # 1e-7 degree offsets; largest difference seen 1.9e-4 km (2026-11-25)
    rng = np.random.default_rng(3)
    lat = rng.uniform(-80, 80, 20000)
    lon = rng.uniform(-180, 180, 20000)
    eps = rng.uniform(-1e-7, 1e-7, (2, 20000))
    got = haversine(lat, lon, -lat + eps[0], lon - 180 + eps[1])
    assert np.all(np.isfinite(got))
    np.testing.assert_allclose(got, np.pi * 6371, rtol=0, atol=1e-3)


def test_haversine_nan_input_stays_nan():
    assert np.isnan(haversine(np.nan, 0, 10, 10))
    assert np.isnan(haversine(0, 0, 10, np.nan))


def test_km2deg_one_degree():
    assert km2deg(6371 * np.pi / 180) == pytest.approx(1, rel=1e-15)
