"""Limb selection and screening: boxes on valid positions and across the
dateline (Q12), extras that would replace the loader's own variables, and
the file's own points screened before interpolation (Q13).

Synthetic Misc files, so every test runs in CI. Known answers: the
reviewer's example (heights 20, 22, 24 km, temperatures 250, 50, 250 K,
interpolated to 21 and 23 km gave 150 K with flag 0) and a finite sentinel
longitude of -999 (which became 81 E after unwrapping and wrapping).
"""

import numpy as np
import pytest
from scipy.io import savemat

from atmosgw.curate.data import states, validate
from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs import load_limb_sounders


def _write(folder, temp, lat=10.0, lon=20.0, heights=(20.0, 22.0, 24.0), **extra):
    alt = np.atleast_2d(np.asarray(heights, dtype=float))
    temp = np.atleast_2d(np.asarray(temp, dtype=float))
    alt = np.broadcast_to(alt, temp.shape).copy()
    content = {"Alt": alt, "Temp": temp,
               "Lat": np.broadcast_to(np.asarray(lat, float), temp.shape).copy(),
               "Lon": np.broadcast_to(np.asarray(lon, float), temp.shape).copy(),
               "Time": np.full(temp.shape, datenum(2007, 1, 1) + 0.5), **extra}
    savemat(folder / "mine_2007d001.mat", content)


def _load(folder, **options):
    return load_limb_sounders("2007-01-01", "Misc", misc_path=folder,
                              misc_identifier="mine", **options)


def test_q13_bad_file_point_not_blended_into_a_good_value(tmp_path):
    _write(tmp_path, [250.0, 50.0, 250.0])
    ds = _load(tmp_path, height_scale=[21.0, 23.0])
    assert np.all(np.isnan(ds.temperature.values))
    assert np.all(states.has_state(ds.temperature_flag.values,
                                   "outside_physical_range"))
    validate(ds, "profile")
    # keep_outliers: the values are kept, the states still set
    ds = _load(tmp_path, height_scale=[21.0, 23.0], keep_outliers=True)
    np.testing.assert_array_equal(ds.temperature.values, [[150.0, 150.0]])
    assert np.all(states.has_state(ds.temperature_flag.values,
                                   "outside_physical_range"))
    # a height that uses only good points is unaffected
    _write(tmp_path, [250.0, 240.0, 230.0, 50.0], heights=(20, 22, 24, 26))
    ds = _load(tmp_path, height_scale=[21.0, 23.0, 25.0])
    np.testing.assert_array_equal(ds.temperature.values[0, :2], [245.0, 235.0])
    assert np.all(ds.temperature_flag.values[0, :2] == 0)
    assert np.isnan(ds.temperature.values[0, 2])


def test_q13_sentinel_longitude_not_wrapped_into_a_position(tmp_path):
    _write(tmp_path, [230.0, 230.0, 230.0], lon=np.array([20.0, -999.0, 20.0]))
    ds = _load(tmp_path, height_scale=[21.0, 23.0])
    assert np.all(np.isnan(ds.longitude.values))
    assert np.all(states.has_state(ds.temperature_flag.values, "bad_geolocation"))


def test_q13_0_to_360_longitudes_are_not_bad(tmp_path):
    _write(tmp_path, [230.0, 230.0, 230.0], lon=np.array([359.0, 359.5, 0.5]))
    ds = _load(tmp_path, height_scale=[21.0, 23.0])
    np.testing.assert_allclose(ds.longitude.values, [[-0.75, 0.0]], atol=1e-12)
    assert np.all(ds.temperature_flag.values == 0)


def test_q12_dateline_and_greenwich_boxes(tmp_path):
    # profile 0 crosses the dateline; profile 1 is at Greenwich
    lon = np.array([[179.0, -179.0, -178.0], [0.5, 0.5, 0.5]])
    _write(tmp_path, np.full((2, 3), 230.0), lon=lon)
    assert _load(tmp_path, original_z=True, lon_range=(-1, 1)).sizes["profile"] == 1
    ds = _load(tmp_path, original_z=True, lon_range=(170, -170))
    assert ds.sizes["profile"] == 1 and ds.longitude.values[0, 0] == 179.0


@pytest.mark.parametrize("name", ["temperature", "Temp", "time", "pressure_flag",
                                  "PresEstimated"])
def test_extras_may_not_replace_the_loaders_variables(tmp_path, name):
    _write(tmp_path, [250.0, 250.0, 250.0], temperature=np.full((1, 3), 299.0))
    with pytest.raises(ValueError, match="would replace"):
        _load(tmp_path, additional_vars=[name])


def test_repeated_extras_loaded_once(tmp_path):
    _write(tmp_path, [250.0, 250.0, 250.0], Extra=np.full((1, 3), 7.0))
    ds = _load(tmp_path, additional_vars=["Extra", "Extra"], original_z=True)
    np.testing.assert_array_equal(ds["Extra"].values, [[7.0, 7.0, 7.0]])
