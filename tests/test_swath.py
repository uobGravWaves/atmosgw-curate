"""Tests for atmosgw.curate.swath.

MATLAB outputs: tests/data/matlab_swath_fixtures.mat (see
tests/matlab/make_swath_fixtures.m): Mapping Toolbox ``distance``, and
``prep_airs_3d`` (corwin365/MatlabFunctions, recorded in the file) regridding
three synthetic swaths with InputStruct, without detrending.

Tolerances against MATLAB (trig functions round differently in the last bits
between MATLAB and NumPy, and MATLAB's ``distance`` may arrange its formula
differently):
- arc and azimuth: 1e-9 degrees absolute (0.1 mm on the ground);
- regridded latitude and longitude: 1e-9 degrees;
- regridded temperature: 1e-9 K; time: 1e-9 days (0.1 ms; datenums of about
  7.3e5 carry about 1e-10 days of rounding);
- grid spacing: 1e-9 km;
- vertical grid: exact.
NaN patterns must match exactly.
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from atmosgw.curate import swath

FIXTURES = Path(__file__).parent / "data" / "matlab_swath_fixtures.mat"
DEG, KELVIN, DAYS, KM = 1e-9, 1e-9, 1e-9, 1e-9
CALLS = {  # MATLAB options -> swath.regrid options (NXT 16, NAT 12 in every call)
    "linear": {},
    "vertical": {},
    "spacing": {"dx_cross": 40, "dx_along": 30},
    "nearest": {"method": "nearest"},
}


@pytest.fixture(scope="module")
def matlab():
    return loadmat(FIXTURES, squeeze_me=True, struct_as_record=False)


def python_prep(inp, call):
    """The prep_airs_3d steps covered here, on a MATLAB InputStruct."""
    temp = np.transpose(inp.ret_temp, (1, 2, 0))  # z x XT x AT -> XT x AT x z
    time = 730486.0 + inp.l1_time / 86400  # datenum(2000,1,1) + seconds
    options = {"n_cross": 16, "n_along": 12}
    if "dx_cross" in CALLS[call]:
        options = {}
    options.update(CALLS[call])
    out = swath.regrid(inp.l1_lat, inp.l1_lon, {"Temp": temp, "Time": time}, **options)
    z = np.asarray(inp.ret_z, dtype=float)
    temp = out["fields"]["Temp"]
    if call == "vertical":
        new_z = swath.regular_levels(z, 3)
        temp = swath.interpolate_levels(z, temp, new_z)
        z = new_z
    return out, temp, z


def assert_close(got, want, atol, what):
    got, want = np.asarray(got, dtype=float), np.asarray(want, dtype=float)
    assert got.shape == want.shape, f"{what}: shape {got.shape} vs {want.shape}"
    np.testing.assert_array_equal(np.isnan(got), np.isnan(want), err_msg=what)
    np.testing.assert_allclose(got, want, rtol=0, atol=atol, equal_nan=True,
                               err_msg=what)


def test_arc_distance_against_matlab(matlab):
    a, b = matlab["DistA"], matlab["DistB"]
    arc, az = swath.arc_distance(a[:, 0], a[:, 1], b[:, 0], b[:, 1])
    assert_close(arc, matlab["DistArc"], DEG, "arc")
    assert_close(az, matlab["DistAz"], DEG, "azimuth")


@pytest.mark.parametrize("geometry", ["mid", "dateline", "turning"])
@pytest.mark.parametrize("call", list(CALLS))
def test_regrid_against_matlab(matlab, geometry, call):
    case = getattr(matlab["Swath"], geometry)
    want = getattr(case, call)
    out, temp, z = python_prep(case.In, call)
    assert_close(out["lat"], want.Lat, DEG, "latitude")
    dlon = np.mod(out["lon"] - want.Lon + 180, 360) - 180
    assert_close(dlon, np.zeros_like(dlon), DEG, "longitude")
    assert_close(out["fields"]["Time"], want.Time, DAYS, "time")
    assert_close(temp, want.Temp, KELVIN, "temperature")
    np.testing.assert_array_equal(z, want.Z)
    assert_close(out["spacing"], np.asarray(want.Spacing)[:2], KM, "spacing")


def test_regrid_longitude_rotation():
    """Shifting every longitude by c shifts the output by c and changes nothing else."""
    lat, lon = synthetic_swath(-10, 120, 15)
    temp = np.random.default_rng(1).normal(250, 5, lat.shape + (3,))
    base = swath.regrid(lat, lon, {"t": temp}, n_cross=20, n_along=30)
    for c in (0.0, 180.0 - float(np.median(lon)), -137.25):
        rot = np.mod(lon + c + 180, 360) - 180
        out = swath.regrid(lat, rot, {"t": temp}, n_cross=20, n_along=30)
        dlon = np.mod(out["lon"] - (base["lon"] + c) + 180, 360) - 180
        assert np.max(np.abs(dlon)) < 1e-9
        assert np.max(np.abs(out["lat"] - base["lat"])) < 1e-9
        np.testing.assert_allclose(out["fields"]["t"], base["fields"]["t"],
                                   rtol=0, atol=1e-9)


def test_regrid_no_jump_at_dateline():
    lat, lon = synthetic_swath(5, 179.0, 80)  # crosses the dateline along track
    out = swath.regrid(lat, lon, {}, n_cross=30, n_along=60)
    m = out["lat"].shape[0] // 2
    step, _ = swath.arc_distance(out["lat"][m, :-1], out["lon"][m, :-1],
                                 out["lat"][m, 1:], out["lon"][m, 1:])
    step = swath._deg2km(step)
    assert np.max(step) / np.median(step) < 1.01
    assert np.all((out["lon"] >= -180) & (out["lon"] <= 180))


def test_regrid_regular_spacing_and_range():
    lat, lon = synthetic_swath(50, -30, 10)
    dist = swath.swath_distances(lat, lon)
    out = swath.regrid(lat, lon, {}, dx_cross=25, dx_along=12)
    np.testing.assert_allclose(np.diff(out["cross"]), 25, atol=1e-9)
    np.testing.assert_allclose(np.diff(out["along"]), 12, atol=1e-9)
    assert out["cross"][0] == dist["cross"].min()
    assert out["cross"][-1] <= dist["cross"].max()
    assert out["spacing"] == pytest.approx((25, 12))


def test_regrid_linear_field_exact_and_nan_spreads():
    """A field linear in distance is reproduced; a NaN reaches only its cells."""
    lat, lon = synthetic_swath(30, 0, 20)
    dist = swath.swath_distances(lat, lon)
    xx, aa = np.meshgrid(dist["cross"], dist["along"],
                         indexing="ij")
    out = swath.regrid(lat, lon, {"f": 2 * xx - 0.5 * aa}, n_cross=13, n_along=17)
    gx, ga = np.meshgrid(out["cross"], out["along"], indexing="ij")
    np.testing.assert_allclose(out["fields"]["f"], 2 * gx - 0.5 * ga, atol=1e-8)
    f = 2 * xx - 0.5 * aa
    f[5, 5] = np.nan
    out = swath.regrid(lat, lon, {"f": f}, n_cross=13, n_along=17)
    near = ((gx > dist["cross"][4]) & (gx < dist["cross"][6])
            & (ga > dist["along"][4]) & (ga < dist["along"][6]))
    assert np.array_equal(np.isnan(out["fields"]["f"]), near)


def test_regular_levels_and_interpolation():
    airs = np.array([0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48,
                     51, 54, 57, 60, 65, 70, 75, 80, 85, 90], dtype=float)
    np.testing.assert_array_equal(swath.regular_levels(airs), np.arange(0, 91, 3.0))
    # as MATLAB: from the lowest level, whatever the middle level (its comment says
    # "from the middle outwards", but min:step:max starts at the bottom)
    z = np.array([1.0, 4.0, 7.5, 9.0, 13.0])
    np.testing.assert_array_equal(swath.regular_levels(z, 2), np.arange(1, 14, 2.0))
    # linear in height: exact; NaN outside the input levels
    values = 200 + 2 * np.broadcast_to(airs, (3, 4, airs.size))
    out = swath.interpolate_levels(airs, values, np.array([-1.0, 1.5, 62.5, 91.0]))
    want = np.broadcast_to(200 + 2 * np.array([1.5, 62.5]), (3, 4, 2))
    np.testing.assert_allclose(out[..., 1:3], want)
    assert np.all(np.isnan(out[..., [0, 3]]))


def test_unwrap_degrees():
    np.testing.assert_array_equal(swath.unwrap_degrees([170, 179, -179, -170, 175]),
                                  [170, 179, 181, 190, 175 + 360 - 360])
    np.testing.assert_array_equal(swath.unwrap_degrees([10, 20]), [10, 20])


def test_regrid_rejects_bad_options():
    lat, lon = synthetic_swath(0, 0, 5)
    with pytest.raises(ValueError):
        swath.regrid(lat, lon, {}, method="natural")
    with pytest.raises(ValueError):
        swath.regrid(lat, lon, {}, extrapolate="none")


# helpers


def synthetic_swath(lat0, lon0, heading, n_cross=31, n_along=40, step=18.0):
    """A curved swath built with reckon: [cross x along] lat and lon."""
    from atmosgw.curate.geodesy import km2deg, reckon

    mid_lat, mid_lon = np.empty(n_along), np.empty(n_along)
    mid_lat[0], mid_lon[0] = lat0, lon0
    head = heading + np.linspace(0, 8, n_along)
    for k in range(1, n_along):
        mid_lat[k], mid_lon[k] = reckon(mid_lat[k - 1], mid_lon[k - 1],
                                        km2deg(step), head[k - 1])
    offsets = np.sign(np.arange(n_cross) - n_cross // 2) * 15.0 * np.abs(
        np.arange(n_cross) - n_cross // 2) ** 1.05
    lat, lon = reckon(np.broadcast_to(mid_lat, (n_cross, n_along)),
                      np.broadcast_to(mid_lon, (n_cross, n_along)),
                      km2deg(np.broadcast_to(offsets[:, None], (n_cross, n_along))),
                      np.broadcast_to(head + 90, (n_cross, n_along)))
    return np.asarray(lat), np.asarray(lon)
