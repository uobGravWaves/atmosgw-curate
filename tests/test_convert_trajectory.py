"""Tests for to_profiles on radiosonde trajectories (data.convert).

Synthetic trajectory Datasets (checked valid with ``validate(...,
"trajectory")``) and synthetic DWD files (``test_new_loaders.write_sondes``);
no real data.

Tolerances: a profile linear in height is reproduced at the grid levels to
1e-10 K (or m s-1, degrees; a + w (b - a) in float64 rounds to a few 1e-14 at
300 K); an exponential pressure profile, interpolated in log(pressure), to a
relative 1e-12 (log and exp round to ~1e-15; linear interpolation in
pressure would be off by ~1e-5 at these sample spacings, so the test tells
them apart); times to 1 us (rounded to the us). Everything else exactly.
"""

import json

import numpy as np
import pytest
import xarray as xr
from test_data_formats import example_swath, example_time_height
from test_new_loaders import write_sondes

from atmosgw.curate.data import FORMATS, states, to_profiles, validate
from atmosgw.curate.obs import load_trajectories

START = np.datetime64("2020-01-02T04:45", "us")


def _flag(value):
    """Loader-style states: missing between a flight's first and last value,
    no_observation outside."""
    flag = states.new_flags(value.shape)
    for k, row in enumerate(value):
        finite = np.flatnonzero(np.isfinite(row))
        inside = np.zeros(row.shape, dtype=bool)
        if finite.size:
            inside[finite[0]:finite[-1] + 1] = True
        states.set_state(flag[k], np.isnan(row) & inside, "missing")
        states.set_state(flag[k], np.isnan(row) & ~inside, "no_observation")
    return flag


def trajectory(alt, temp, pres="exp", lat=None, lon=None, instrument="DWDsondes"):
    """A valid trajectory Dataset; alt, temp [flight x sample] (km, K); pres
    "exp" (1000 exp(-z / 7) hPa), an array, or None (absent); winds
    u = 2 z, v = -z (m s-1); samples 2 s apart."""
    alt = np.atleast_2d(np.asarray(alt, dtype=float))
    temp = np.atleast_2d(np.asarray(temp, dtype=float))
    n_f, n_s = alt.shape
    ts = ("flight", "sample")
    lat = np.full(alt.shape, 50.0) if lat is None else np.atleast_2d(lat).astype(float)
    lon = np.full(alt.shape, 10.0) if lon is None else np.atleast_2d(lon).astype(float)
    values = {"temperature": (temp, "K"), "eastward_wind": (2 * alt, "m s-1"),
              "northward_wind": (-alt, "m s-1")}
    if isinstance(pres, str):
        values["pressure"] = (1000 * np.exp(-alt / 7), "hPa")
    elif pres is not None:
        values["pressure"] = (np.atleast_2d(pres).astype(float), "hPa")
    variables = {}
    for name, (value, units) in values.items():
        variables[name] = (ts, value, {"units": units, "origin": "measured",
                                       "ancillary_variables": f"{name}_flag"})
        variables[f"{name}_flag"] = (ts, _flag(value), states.flag_attrs())
    variables.update({
        "flight_id": (("flight",), np.array([f"01234-{k}" for k in range(n_f)])),
        "platform_id": (("flight",), np.array(["01234"] * n_f)),
        "platform_name": (("flight",), np.array(["Teststadt"] * n_f)),
        "platform_type": (("flight",), np.array(["Radiosonde TEST"] * n_f)),
    })
    coords = {
        "latitude": (ts, lat, {"units": "degrees_north"}),
        "longitude": (ts, lon, {"units": "degrees_east"}),
        "time": (ts, START + np.broadcast_to(np.arange(n_s) * 2, alt.shape)
                 * np.timedelta64(1, "s")),
        "altitude": (ts, alt, {"units": "km", "comment": "geopotential height"}),
        "start_time": (("flight",), np.full(n_f, START)),
    }
    absent = {"upward_wind": "not_measured", "source_file": "not_requested",
              "original_files": "not_requested"}
    if pres is None:
        absent["pressure"] = "not_in_files"
    absent.update({f"{n}_uncertainty": "not_provided" for n in values})
    attrs = {a: "" for a in FORMATS["trajectory"]["attrs"]}
    attrs.update(format="trajectory", format_version=1, instrument=instrument,
                 product="sounding", characterisation=f"{instrument}-SOUNDING-"
                 "placeholder/0", availability=states.format_availability(absent))
    ds = xr.Dataset(variables, coords=coords, attrs=attrs)
    validate(ds, "trajectory")
    return ds


def ascent(top=3.0, bottom=0.1, seed=0):
    """Irregular, strictly increasing sample heights (km), 3-8 m apart."""
    steps = np.random.default_rng(seed).uniform(0.003, 0.008, 2000)
    z = bottom + np.r_[0, np.cumsum(steps)]
    return z[z <= top]


def lapse(z):
    return 290.0 - 6.5 * z


def test_linear_profile_reproduced_on_grid():
    z = ascent()
    p = to_profiles(trajectory(z, lapse(z)))
    validate(p, "profile")
    k = np.arange(int(np.ceil(z[0] / 0.05)), int(np.floor(z[-1] / 0.05)) + 1)
    np.testing.assert_array_equal(p.altitude.values[0], k * 0.05)  # k * dz exactly
    assert np.all(np.diff(p.altitude.values[0]) > 0)
    g = p.altitude.values[0]
    np.testing.assert_allclose(p.temperature.values[0], lapse(g), rtol=0, atol=1e-10)
    np.testing.assert_allclose(p.eastward_wind.values[0], 2 * g, rtol=0, atol=1e-10)
    assert p.temperature.dtype == np.float64 and p.latitude.dtype == np.float64
    assert np.all(p.temperature_flag.values == 0)
    assert "linearly interpolated" in p.temperature.attrs["comment"]
    assert "log(pressure)" in p.pressure.attrs["comment"]
    assert p.attrs["converted_from"] == "trajectory"
    assert p.attrs["format"] == "profile"
    assert p.attrs["format_version"] == FORMATS["profile"]["version"]
    conversion = json.loads(p.attrs["conversion"])
    assert conversion["dz_km"] == 0.05 and conversion["height_range"] is None
    assert list(p.flight_id.values) == ["01234-0"]
    assert p.start_time.values[0] == START
    # time is linear in the sample index; at each level, between its samples
    sec = np.interp(g, z, np.arange(z.size) * 2.0)
    expected = START + np.round(sec * 1e6).astype(np.int64).astype("timedelta64[us]")
    assert np.all(np.abs(p.time.values[0] - expected) <= np.timedelta64(1, "us"))


def test_float32_input_keeps_its_type():
    z = ascent(top=1.0)
    ds = trajectory(z, lapse(z))
    ds["temperature"] = ds.temperature.astype(np.float32)
    p = to_profiles(ds)
    validate(p, "profile")
    assert p.temperature.dtype == np.float32 and p.pressure.dtype == np.float64


def test_descent_after_burst_is_dropped():
    up, down = ascent(top=5.0), np.linspace(4.99, 2.0, 300)
    z = np.r_[up, down]
    t = np.r_[lapse(up), np.full(down.size, 250.0)]  # (descent: other values)
    p = to_profiles(trajectory(z, t))
    g = p.altitude.values[0]
    assert g[-1] == np.floor(up[-1] / 0.05) * 0.05
    np.testing.assert_allclose(p.temperature.values[0], lapse(g), rtol=0, atol=1e-10)


def test_wiggles_use_the_running_maximum():
    z = np.array([1.00, 1.11, 1.04, 1.08, 1.13, 1.12, 1.21, 1.32])
    t = np.array([220.0, 222.0, 300.0, 300.0, 224.0, 300.0, 226.0, 228.0])
    p = to_profiles(trajectory(z, t))
    g = p.altitude.values[0]
    np.testing.assert_allclose(g, [1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3])
    kept_z, kept_t = z[[0, 1, 4, 6, 7]], t[[0, 1, 4, 6, 7]]
    np.testing.assert_allclose(p.temperature.values[0], np.interp(g, kept_z, kept_t),
                               rtol=0, atol=1e-10)
    assert np.all(p.temperature.values[0] < 230)  # (no 300 K sample used)


def test_no_extrapolation_outside_each_flight():
    za, zb = ascent(top=1.0, bottom=0.2), ascent(top=2.0, bottom=0.5, seed=1)
    n = max(za.size, zb.size)
    alt = np.full((2, n), np.nan)
    alt[0, :za.size], alt[1, :zb.size] = za, zb
    temp = lapse(alt)
    ds = trajectory(alt, temp)
    ds["time"] = ds.time.where(np.isfinite(alt))  # (padding: NaT)
    validate(ds, "trajectory")
    p = to_profiles(ds)
    validate(p, "profile")
    g = p.altitude.values[0]
    assert g[0] == 0.2 and g[-1] == np.floor(zb[-1] / 0.05) * 0.05
    for k, z in enumerate((za, zb)):
        outside = (g < z[0]) | (g > z[-1])
        assert outside.any()
        for name in ("temperature", "pressure", "eastward_wind"):
            assert np.all(np.isnan(p[name].values[k, outside]))
            assert np.all(np.isfinite(p[name].values[k, ~outside]))
            flag = p[f"{name}_flag"].values[k]
            assert np.all(states.has_state(flag[outside], "no_observation"))
            assert not states.has_state(flag, "missing").any()
        assert np.all(np.isnat(p.time.values[k, outside]))
        assert np.all(np.isnan(p.latitude.values[k, outside]))
        assert np.all(np.isnan(p.longitude.values[k, outside]))
    assert "no_observation" in p.attrs["state_counts"]


def test_a_missing_sample_is_not_bridged():
    z = np.array([0.10, 0.13, 0.17, 0.22, 0.26, 0.31, 0.36])
    t = lapse(z)
    t[3] = np.nan  # (0.22 km)
    p = to_profiles(trajectory(z, t))
    g = p.altitude.values[0]
    gap = (g > z[2]) & (g < z[4])
    assert gap.sum() == 2  # (0.20, 0.25)
    assert np.all(np.isnan(p.temperature.values[0, gap]))
    assert np.all(states.has_state(p.temperature_flag.values[0, gap], "missing"))
    np.testing.assert_allclose(p.temperature.values[0, ~gap], lapse(g[~gap]),
                               rtol=0, atol=1e-10)
    assert np.all(np.isfinite(p.pressure.values[0]))  # (other variables kept)


def test_removal_state_of_a_sample_is_carried(tmp_path):
    # write_sondes: on 2020-01-02 the 05 ascent's sample 5 (0.15 km) is 50 K,
    # removed as outside_physical_range; sample 7 has latitude 0.00 (a fill)
    write_sondes(tmp_path, n_days=2)
    ds = load_trajectories("2020-01-02", "DWDsondes", file_path=tmp_path,
                           file_source=True)
    p = to_profiles(ds)
    validate(p, "profile")
    g = p.altitude.values[0]
    np.testing.assert_allclose(g, np.arange(2, 14) * 0.05)
    at = int(np.flatnonzero(np.isclose(g, 0.15))[0])
    assert np.isnan(p.temperature.values[0, at])
    flag = p.temperature_flag.values[0, at]
    assert states.has_state(flag, "outside_physical_range")
    assert not states.has_state(flag, "missing")
    # elsewhere 10 - k degrees C at 0.1 + 0.01 k km: linear in height
    other = np.arange(g.size) != at
    np.testing.assert_allclose(p.temperature.values[0, other],
                               283.15 - (g[other] - 0.1) / 0.01, rtol=0, atol=1e-10)
    assert np.isnan(p.temperature.values[1, at])  # (the 11 UTC ascent too)
    # per profile, provenance and availability
    assert list(p.flight_id.values) == list(ds.flight_id.values)
    assert list(p.reference_time.values) == list(ds.reference_time.values)
    assert p.source_file.dims == ("profile", "level")
    assert np.all(p.source_file.values == ds.source_file.values[:, None])
    assert list(p.original_files.values) == list(ds.original_files.values)
    avail = states.parse_availability(p.attrs["availability"])
    assert avail["source_profile"] == "not_provided"
    assert avail["upward_wind"] == "not_measured"
    assert avail["temperature_uncertainty"] == "not_provided"
    assert p.pressure.attrs["origin"] == "measured"


def test_dateline_longitude():
    z = np.linspace(0.1, 1.0, 91)  # (10 m apart)
    lon = 179.95 + 0.1 * (z - 0.1) / 0.9 * 9  # (179.95 to 180.85 east)
    p = to_profiles(trajectory(z, lapse(z), lon=np.where(lon > 180, lon - 360, lon)))
    validate(p, "profile")
    g = p.altitude.values[0]
    expected = 179.95 + 0.1 * (g - 0.1) / 0.9 * 9
    expected = np.where(expected > 180, expected - 360, expected)
    out = p.longitude.values[0]
    assert np.all((out >= -180) & (out <= 180))
    np.testing.assert_allclose(out, expected, rtol=0, atol=1e-9)
    assert (out > 179).any() and (out < -179).any()


def test_pressure_interpolated_in_log_pressure():
    z = np.arange(0.1, 20.0, 0.137)  # (137 m apart: linear would err by ~1e-5)
    p = to_profiles(trajectory(z, lapse(z)))
    g = p.altitude.values[0]
    np.testing.assert_allclose(p.pressure.values[0], 1000 * np.exp(-g / 7),
                               rtol=1e-12, atol=0)
    linear = np.interp(g, z, 1000 * np.exp(-z / 7))
    assert np.max(np.abs(linear / p.pressure.values[0] - 1)) > 1e-6
    assert np.all(p.pressure_flag.values == 0)


def test_absent_pressure_is_estimated():
    from atmosgw.curate.atmosphere import height_to_pressure

    z = ascent(top=1.0)
    p = to_profiles(trajectory(z, lapse(z), pres=None))
    validate(p, "profile")
    np.testing.assert_array_equal(p.pressure.values,
                                  height_to_pressure(p.altitude.values))
    assert p.pressure.attrs["origin"] == "derived"
    assert states.has_state(p.pressure_flag.values, "derived").all()


def test_height_range_and_dz():
    z = ascent(top=2.0, bottom=0.12)
    ds = trajectory(z, lapse(z))
    p = to_profiles(ds, dz_km=0.1, height_range=(0, 3))
    validate(p, "profile")
    np.testing.assert_array_equal(p.altitude.values[0], np.arange(31) * 0.1)
    g = p.altitude.values[0]
    covered = (g >= z[0]) & (g <= z[-1])
    np.testing.assert_allclose(p.temperature.values[0, covered], lapse(g[covered]),
                               rtol=0, atol=1e-10)
    assert np.all(states.has_state(p.temperature_flag.values[0, ~covered],
                                   "no_observation"))
    conversion = json.loads(p.attrs["conversion"])
    assert conversion["dz_km"] == 0.1 and conversion["height_range"] == [0.0, 3.0]
    q = to_profiles(ds, dz_km=0.25)
    assert 1.75 < z[-1] < 2.0
    np.testing.assert_array_equal(q.altitude.values[0], np.arange(1, 8) * 0.25)
    for bad in ({"dz_km": 0}, {"dz_km": -0.05}, {"height_range": (3, 1)},
                {"height_range": 5}):
        with pytest.raises(ValueError):
            to_profiles(ds, **bad)


def test_aircraft_refused():
    z = ascent(top=1.0)
    for name in ("IAGOS", "FAAM"):
        with pytest.raises(ValueError, match="only radiosonde"):
            to_profiles(trajectory(z, lapse(z), instrument=name))


def test_options_refused_for_swath_and_time_height():
    for ds in (example_swath(), example_time_height()):
        validate(to_profiles(ds), "profile")  # (default call unchanged)
        with pytest.raises(TypeError, match="dz_km"):
            to_profiles(ds, dz_km=0.1)
        with pytest.raises(TypeError, match="height_range"):
            to_profiles(ds, height_range=(0, 10))
