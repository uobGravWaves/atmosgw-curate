"""Known-answer tests of the SABER, HIRDLS, MIPAS and SOFIE readers and the
AIRS daily-file layout, on small files written here in each product's layout,
so that they run in CI (the real-data tests run only where the data are).

Each test checks what a reader could get wrong without failing: the time
epoch and units, the position convention (longitude wrap, the HIRDLS
along-track shift and its sign), height and pressure units, and which events
belong to the day. Expected values are worked out by hand from the file
contents, not with the package's own functions (except height_to_pressure,
to write a plausible pressure).
"""

import h5py
import netCDF4
import numpy as np
import pytest
from test_nadir_known import seconds, write_granule

from atmosgw.curate.atmosphere import height_to_pressure
from atmosgw.curate.data import states, validate
from atmosgw.curate.obs import load_limb_sounders, load_nadir_sounders

KM_PER_DEG = 6371.0 * np.pi / 180  # (the package's sphere)


def _times(ds):
    return ds.time.values[:, 0]


# ------------------------------------------------------------------ SABER


def _saber(tmp_path):
    """January 2007 file, 3 events x 5 levels (GATS Level 2A layout)."""
    folder = tmp_path / "2007"
    folder.mkdir()
    alt = np.array([20.0, 25.0, 30.0, 35.0, 40.0])
    date = np.array([2007001, 2007002, 2007002], dtype=np.int32)
    # ms since midnight: 13:53:20, 01:00, 12:00
    ms = np.array([50_000_000, 3_600_000, 43_200_000], dtype=float)
    path = folder / "SABER_Temp_O3_H2O_January2007_v2.0.nc"
    with netCDF4.Dataset(path, "w") as f:
        f.createDimension("event", 3)
        f.createDimension("altitude", 5)
        f.Data_Product_Version = "02.07"

        def var(name, value, dims=("event", "altitude"), kind="f4"):
            v = f.createVariable(name, kind, dims)
            v.missing_value = np.array(-999, dtype=kind)
            v[:] = value
        var("date", date, ("event",), "i4")
        var("time", np.tile(ms[:, None], (1, 5)), kind="i4")
        tpalt = np.tile(alt, (3, 1))
        tpalt[2, 4] = 2000.0  # (a bad height: above 1000 km becomes NaN)
        var("tpaltitude", tpalt)
        var("tplatitude", np.full((3, 5), 10.0))
        var("tplongitude", np.tile(np.array([[10.0], [20.0], [200.0]]), (1, 5)))
        var("pressure", np.tile(height_to_pressure(alt), (3, 1)))
        var("ktemp", 200.0 + np.arange(15.0).reshape(3, 5))
    return alt


def test_saber_events_times_positions(tmp_path):
    alt = _saber(tmp_path)
    ds = load_limb_sounders("2007-01-02", "SABER", file_path=tmp_path,
                            original_z=True, keep_outliers=True, file_source=True)
    # events 2 and 3 are dated 2 January: 01:00 and 12:00 (date + ms since midnight)
    assert ds.source_profile.values[:, 0].tolist() == [2, 3]
    np.testing.assert_array_equal(
        _times(ds), np.array(["2007-01-02T01:00", "2007-01-02T12:00"],
                             dtype="datetime64[us]"))
    np.testing.assert_array_equal(ds.longitude.values[:, 0], [20.0, -160.0])  # (0-360)
    np.testing.assert_array_equal(ds.altitude.values[0], alt)
    assert np.isnan(ds.altitude.values[1, 4])  # (2000 km)
    np.testing.assert_array_equal(ds.temperature.values[0], 205.0 + np.arange(5.0))
    np.testing.assert_allclose(ds.pressure.values[0], height_to_pressure(alt),
                               rtol=1e-6)  # (float32 in the file)
    assert ds.attrs["source_versions"] == "Data_Product_Version 02.07 (1 file)"


# ----------------------------------------------------------------- HIRDLS


def _hirdls(tmp_path, scan_up=(1, 0, 1), precision=1.0):
    """3 profiles going north along 0 E, 5 levels (HDF-EOS5 layout)."""
    folder = tmp_path / "2007"
    folder.mkdir(exist_ok=True)
    path = folder / "HIRDLS-Aura_L2_v07-00-20-c01_2007d001.he5"
    alt_km = np.array([20.0, 25.0, 30.0, 35.0, 40.0])
    tai = (np.datetime64("2007-01-01T06:00") - np.datetime64("1993-01-01")) \
        / np.timedelta64(1, "s") + np.array([0.0, 60.0, 120.0])
    with h5py.File(path, "w") as f:
        geo = "HDFEOS/SWATHS/HIRDLS/Geolocation Fields"
        dat = "HDFEOS/SWATHS/HIRDLS/Data Fields"
        f[f"{geo}/Latitude"] = np.array([0.0, 1.0, 2.0], np.float32)
        f[f"{geo}/Longitude"] = np.zeros(3, np.float32)
        f[f"{geo}/Altitude"] = np.tile(alt_km * 1000, (3, 1)).astype(np.float32)
        f[f"{geo}/Pressure"] = height_to_pressure(alt_km).astype(np.float32)
        f[f"{geo}/Time"] = tai
        f[f"{geo}/ScanUpFlag"] = np.array(scan_up, np.float32)
        f[f"{dat}/Temperature"] = (220.0 + np.arange(15.0).reshape(3, 5)).astype(
            np.float32)
        prec = np.full((3, 5), precision, np.float32)
        prec[0, 0] = -1.0  # (mostly a priori)
        f[f"{dat}/TemperaturePrecision"] = prec
    return alt_km


def test_hirdls_time_heights_and_geolocation_shift(tmp_path):
    alt = _hirdls(tmp_path)
    ds = load_limb_sounders("2007-01-01", "HIRDLS", file_path=tmp_path,
                            original_z=True)
    np.testing.assert_array_equal(
        _times(ds), np.array(["2007-01-01T06:00", "2007-01-01T06:01",
                              "2007-01-01T06:02"], dtype="datetime64[us]"))
    np.testing.assert_array_equal(ds.altitude.values[0], alt)  # (m -> km)
    # each point moved 0.6 km per km from 30 km along the track (north here):
    # forwards for up-scans, backwards for down-scans [R2]
    shift = 0.6 * (alt - 30) / KM_PER_DEG
    sign = np.array([1, -1, 1])[:, None]
    expected = np.array([0.0, 1.0, 2.0])[:, None] + sign * shift
    np.testing.assert_allclose(ds.latitude.values, expected, atol=1e-9)
    np.testing.assert_allclose(ds.longitude.values, 0.0, atol=1e-9)
    # negative precision: kept, marked a_priori_dominated (Q17)
    flag = ds.temperature_flag.values
    assert states.has_state(flag[0, 0], "a_priori_dominated")
    assert ds.temperature[0, 0] == 220
    assert not states.has_state(flag, "a_priori_dominated")[0, 1:].any()
    validate(ds, "profile")


# ------------------------------------------------------------------ MIPAS


def test_mipas_time_from_hours_and_values(tmp_path):
    folder = tmp_path / "2007"
    folder.mkdir()
    with netCDF4.Dataset(folder / "mipas_envisat_20070102_nom_oxf_v1.40.nc", "w") as f:
        f.createDimension("profiles", 2)
        f.createDimension("levels", 4)
        f.history = "L2 processing version 1.40"
        for name, value in (("TEM", 220.0 + np.arange(8.0).reshape(2, 4)),
                            ("PRE", np.tile(height_to_pressure([20, 30, 40, 50]),
                                            (2, 1))),
                            ("hgt", np.tile([20.0, 30.0, 40.0, 50.0], (2, 1)))):
            v = f.createVariable(name, "f4", ("profiles", "levels"), fill_value=-999.0)
            v[:] = value
        for name, value in (("latitude", [-45.0, 60.0]), ("longitude", [170.0, -30.0]),
                            ("time", [0.5, 23.25])):  # hours since the file's day
            f.createVariable(name, "f4", ("profiles",))[:] = value
        f["TEM"][1, 3] = -999.0  # (a fill value)
    ds = load_limb_sounders("2007-01-02", "MIPAS", file_path=tmp_path, original_z=True)
    np.testing.assert_array_equal(
        _times(ds), np.array(["2007-01-02T00:30", "2007-01-02T23:15"],
                             dtype="datetime64[us]"))
    np.testing.assert_array_equal(ds.latitude.values[:, 0], [-45.0, 60.0])
    np.testing.assert_array_equal(ds.altitude.values[1], [20.0, 30.0, 40.0, 50.0])
    assert np.isnan(ds.temperature.values[1, 3])
    assert states.has_state(ds.temperature_flag.values[1, 3], "no_observation")
    assert ds.attrs["source_versions"] == "L2 processing version 1.40 (1 file)"


# ------------------------------------------------------------------ SOFIE


def test_sofie_time_from_hours_and_endorsed_heights(tmp_path):
    folder = tmp_path / "2007"
    folder.mkdir()
    alt = np.array([90.0, 50.0, 20.0, 10.0])  # (the files' order: top down)
    with netCDF4.Dataset(folder / "SOFIE_Level2_2007152_01.3.nc", "w") as f:
        f.createDimension("event", 2)
        f.createDimension("altitude", 4)
        f.DP_Version = "01.3"
        temp = np.tile([180.0, 260.0, 220.0, 230.0], (2, 1))
        for name, value in (("Temperature", temp),
                            ("Pressure", np.tile(height_to_pressure(alt), (2, 1)))):
            v = f.createVariable(name, "f4", ("event", "altitude"))
            v.missing_value = np.float32(-1e24)
            v[:] = value
        f.createVariable("Altitude", "f4", ("altitude",))[:] = alt
        f.createVariable("Latitude_83km", "f4", ("event",))[:] = [70.0, -70.0]
        f.createVariable("Longitude_83km", "f4", ("event",))[:] = [10.0, 190.0]
        v = f.createVariable("Time_UT", "f4", ("event",))  # hours, though "msec"
        v.units = "msec"
        v[:] = [6.5, 18.0]
    ds = load_limb_sounders("2007-06-01", "SOFIE", file_path=tmp_path,
                            original_z=True, keep_outliers=True)
    np.testing.assert_array_equal(
        _times(ds), np.array(["2007-06-01T06:30", "2007-06-01T18:00"],
                             dtype="datetime64[us]"))
    np.testing.assert_array_equal(ds.longitude.values[:, 0], [10.0, -170.0])
    np.testing.assert_array_equal(ds.altitude.values[0], alt)
    # outside the endorsed 15-88 km: instrument_flag (Q17)
    bad = states.has_state(ds.temperature_flag.values[0], "instrument_flag")
    np.testing.assert_array_equal(bad, [True, False, False, True])


# -------------------------------------------------------- AIRS daily files


def test_airs_daily_file_equals_granule_files(tmp_path):
    # the same three granules as granule files and as one daily file (the
    # group's layout: dimension "time" over granules) give the same swath
    granules = tmp_path / "granules"
    daily = tmp_path / "daily"
    for g in (1, 2, 3):
        write_granule(granules, g, seconds(g))
    arrays = {}
    for g in (1, 2, 3):
        name_g = granules / "2007" / "001" / f"airs_2007_001_{g:03d}.nc"
        with netCDF4.Dataset(name_g) as f:
            for name in ("l1_lat", "l1_lon", "l1_time", "ret_temp", "ret_z"):
                arrays.setdefault(name, []).append(np.asarray(f[name][:]))
    (daily / "2007").mkdir(parents=True)
    with netCDF4.Dataset(daily / "2007" / "airs_2007d001.nc", "w") as f:
        f.createDimension("time", 3)
        f.createDimension("ALONG", 135)
        f.createDimension("ACROSS", 90)
        f.createDimension("ALTITUDE", 27)
        for name in ("l1_lat", "l1_lon", "l1_time"):
            f.createVariable(name, "f8", ("time", "ALONG", "ACROSS"))[:] = np.stack(
                arrays[name])
        dims = ("time", "ALONG", "ACROSS", "ALTITUDE")
        f.createVariable("ret_temp", "f4", dims)[:] = (
            np.stack(arrays["ret_temp"]))
        f.createVariable("ret_z", "f4", ("ALTITUDE",))[:] = arrays["ret_z"][0]
    request = [("2007-01-01", g) for g in (1, 2, 3)]
    a = load_nadir_sounders(None, "AIRS", granules=request, file_path=granules)
    b = load_nadir_sounders(None, "AIRS", granules=request, file_path=daily)
    for name in ("temperature", "latitude", "longitude", "time", "granule",
                 "distance_along"):
        np.testing.assert_array_equal(a[name].values, b[name].values, err_msg=name)
    assert a.attrs["data"] == "raw" and b.attrs["data"] == "preprocessed"
    assert b.attrs["files_read"] == "airs_2007d001.nc"


@pytest.mark.parametrize("granule", [2, 3])
def test_airs_daily_file_granule_rows(tmp_path, granule):
    # granule g is row g - 1 of the daily file (a 0- or 1-based slip gives
    # the neighbour's times)
    granules, daily = tmp_path / "granules", tmp_path / "daily"
    for g in (1, 2, 3):
        write_granule(granules, g, seconds(g))
    (daily / "2007").mkdir(parents=True)
    with netCDF4.Dataset(daily / "2007" / "airs_2007d001.nc", "w") as f:
        f.createDimension("time", 3)
        f.createDimension("ALONG", 135)
        f.createDimension("ACROSS", 90)
        f.createDimension("ALTITUDE", 27)
        stacks = {}
        for g in (1, 2, 3):
            with netCDF4.Dataset(granules / "2007" / "001" /
                                 f"airs_2007_001_{g:03d}.nc") as src:
                for name in ("l1_lat", "l1_lon", "l1_time", "ret_temp"):
                    stacks.setdefault(name, []).append(np.asarray(src[name][:]))
                z = np.asarray(src["ret_z"][:])
        for name in ("l1_lat", "l1_lon", "l1_time"):
            f.createVariable(name, "f8", ("time", "ALONG", "ACROSS"))[:] = np.stack(
                stacks[name])
        dims = ("time", "ALONG", "ACROSS", "ALTITUDE")
        f.createVariable("ret_temp", "f4", dims)[:] = (
            np.stack(stacks["ret_temp"]))
        f.createVariable("ret_z", "f4", ("ALTITUDE",))[:] = z
    ds = load_nadir_sounders(None, "AIRS", granules=[("2007-01-01", granule)],
                             file_path=daily, original_grid=True)
    first = np.datetime64("2000-01-01", "us") + np.timedelta64(
        int(round(seconds(granule) * 1e6)), "us")
    assert abs(ds.time.values[0, 0] - first) < np.timedelta64(1, "ms")
