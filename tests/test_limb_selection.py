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


@pytest.mark.parametrize("method", ["pchip", "spline", "makima"])
def test_q13_smooth_methods_with_bad_points_do_not_crash(tmp_path, method):
    _write(tmp_path, [250.0, 50.0, 50.0, 50.0], heights=(20, 22, 24, 26))
    ds = _load(tmp_path, height_scale=[21.0, 25.0], interp_type=method)
    assert np.all(np.isnan(ds.temperature.values))
    # a bad point at the top: the heights its extrapolation reaches are marked
    _write(tmp_path, [250.0, 240.0, 230.0, 50.0], heights=(20, 22, 24, 26))
    ds = _load(tmp_path, height_scale=[21.0, 25.0, 27.0], interp_type=method)
    assert np.isfinite(ds.temperature.values[0, 0])
    assert np.all(np.isnan(ds.temperature.values[0, 1:]))
    assert np.all(states.has_state(ds.temperature_flag.values[0, 1:],
                                   "outside_physical_range"))


@pytest.mark.parametrize("method", ["next", "previous"])
def test_q13_whole_number_extras_follow_the_temperature(tmp_path, method):
    _write(tmp_path, [250.0, 240.0, 50.0, 230.0, 220.0], heights=(20, 22, 24, 26, 28),
           Flag=np.full((1, 5), 1.0))
    ds = _load(tmp_path, height_scale=np.arange(20.5, 28, 0.5), interp_type=method,
               additional_vars=["Flag"])
    t, f = ds.temperature.values, ds["Flag"].values
    np.testing.assert_array_equal(np.isnan(t), np.isnan(f))


def _mls(folder, status=(0, 0, 0), quality=(1.5, 1.5, 1.5), convergence=(1.0, 1.0, 1.0),
         precision=None, version="v6.04"):
    """Three MLS profiles on 7 pressure levels (1000 ... 0.0001 hPa)."""
    import h5py
    # (261.016 and 4.6416e-4: the levels as stored in real files, float32)
    pres = np.array([1000, 316, 100, 82.5, 10, 0.001, 0.0001], dtype=np.float32)
    (folder / "2026").mkdir(exist_ok=True)
    path = folder / "2026" / "MLS-Aura_L2GP-Temperature_v06-04-c01_2026d001.he5"
    group = "HDFEOS/SWATHS/Temperature"
    with h5py.File(path, "w") as f:
        f[f"{group}/Data Fields/L2gpValue"] = np.full((3, 7), 230, np.float32)
        f[f"{group}/Data Fields/L2gpPrecision"] = (
            np.full((3, 7), 1.0, np.float32) if precision is None else precision)
        f[f"{group}/Data Fields/Status"] = np.array(status, np.int32)
        f[f"{group}/Data Fields/Quality"] = np.array(quality, np.float32)
        f[f"{group}/Data Fields/Convergence"] = np.array(convergence, np.float32)
        geo = f"{group}/Geolocation Fields"
        f[f"{geo}/Latitude"] = np.array([-10, 0, 10], np.float32)
        f[f"{geo}/Longitude"] = np.array([0, 1, 2], np.float32)
        f[f"{geo}/Pressure"] = pres
        start = ((np.datetime64("2026-01-01") - np.datetime64("1993-01-01"))
                 / np.timedelta64(1, "s"))
        f[f"{geo}/Time"] = start + np.array([0.0, 3600, 7200])
        attrs = f.create_group("HDFEOS/ADDITIONAL/FILE_ATTRIBUTES").attrs
        attrs["PGEVersion"] = np.bytes_(version)
    return pres


def _mls_flags(folder, **options):
    ds = load_limb_sounders("2026-01-01", "MLS", file_path=folder, original_z=True,
                            **options)
    return ds, states.has_state(ds.temperature_flag.values, "instrument_flag")


def test_q15_mls_screening_rules(tmp_path):
    # profile 0 good; profile 1 odd Status; profile 2 Quality 0.5 (good at
    # 83 hPa and above, bad at 100 hPa and below); precision <= 0 at one point
    precision = np.ones((3, 7), np.float32)
    precision[0, 4] = -1.0
    _mls(tmp_path, status=(0, 1, 0), quality=(1.5, 1.5, 0.5), precision=precision)
    ds, bad = _mls_flags(tmp_path)
    # pressure range: 1000, 316 and 0.0001 hPa are outside 261-0.00046 hPa
    np.testing.assert_array_equal(bad[0], [1, 1, 0, 0, 1, 0, 1])
    assert bad[1].all()
    np.testing.assert_array_equal(bad[2], [1, 1, 1, 0, 0, 0, 1])
    assert np.array_equal(np.isnan(ds.temperature.values), bad)
    assert "data quality document" in ds.temperature.attrs["comment"]
    # Convergence must be below 1.03
    _mls(tmp_path, convergence=(1.0, 1.04, 1.0))
    assert _mls_flags(tmp_path)[1][1].all()
    # keep_outliers keeps the values; instrument_qc=False screens nothing
    _mls(tmp_path, status=(0, 1, 0))
    ds, bad = _mls_flags(tmp_path, keep_outliers=True)
    assert bad[1].all() and np.all(ds.temperature.values == 230)
    assert not _mls_flags(tmp_path, instrument_qc=False)[1].any()


def test_q15_mls_unknown_version_not_screened(tmp_path):
    _mls(tmp_path, status=(1, 1, 1), version="v4.23")
    ds, bad = _mls_flags(tmp_path)
    assert not bad.any() and "no MLS screening rules" in ds.temperature.attrs["comment"]


def test_q15_mls_range_ends_as_stored_are_usable(tmp_path):
    _mls(tmp_path)
    import h5py
    path = next((tmp_path / "2026").iterdir())
    with h5py.File(path, "a") as f:
        geo = "HDFEOS/SWATHS/Temperature/Geolocation Fields/Pressure"
        f[geo][...] = np.array([1000, 261.01572, 100, 82.5, 10, 4.6416e-4, 1e-4],
                               np.float32)
    bad = _mls_flags(tmp_path)[1]
    np.testing.assert_array_equal(bad[0], [1, 0, 0, 0, 0, 0, 1])


def _real(instrument, day, **options):
    import os
    import warnings

    from atmosgw.curate.paths import DataDirNotSetError, dataset_dir
    try:
        if not os.path.isdir(dataset_dir(instrument)):
            pytest.skip(f"no {instrument} data")
    except (DataDirNotSetError, FileNotFoundError) as error:
        pytest.skip(f"no {instrument} data: {error}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_limb_sounders(day, instrument, **options)


def test_q17_hirdls_a_priori_marked_not_removed():
    # (real data: about 14 % of HIRDLS values on 2007-01-01 have negative
    # precision, almost all below 20 km)
    ds = _real("HIRDLS", "2007-01-01", original_z=True, keep_outliers=True)
    off = _real("HIRDLS", "2007-01-01", original_z=True, keep_outliers=True,
                instrument_qc=False)
    a_priori = states.has_state(ds.temperature_flag.values, "a_priori_dominated")
    assert 0.05 < a_priori.mean() < 0.3
    np.testing.assert_array_equal(ds.temperature.values, off.temperature.values)
    assert not states.has_state(off.temperature_flag.values,
                                "a_priori_dominated").any()


def test_q17_sofie_outside_endorsed_heights_removed():
    kept = _real("SOFIE", "2007-06-01", original_z=True, keep_outliers=True)
    alt = kept.altitude.values
    bad = states.has_state(kept.temperature_flag.values, "instrument_flag")
    outside = (alt < 15) | (alt > 88)
    assert outside.any()
    np.testing.assert_array_equal(bad, outside)
    ds = _real("SOFIE", "2007-06-01", original_z=True)
    assert np.all(np.isnan(ds.temperature.values[outside]))
    assert np.all(np.isfinite(ds.temperature.values[~outside])
                  == np.isfinite(kept.temperature.values[~outside]))
