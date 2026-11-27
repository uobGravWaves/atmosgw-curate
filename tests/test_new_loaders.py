"""Known-answer tests for the loaders added on 2026-11-27: LIMS, Aeolus and
ROMEX (load_limb_sounders), netCDF sampler output (SWANS) and
``sampled_from`` (limb and nadir), DWD radiosondes (load_trajectories, the
trajectory format) and profile format version 2.

Synthetic files are written to temporary folders, so no real data are
needed; real-data checks are skipped if the files are absent.

Tolerances: values converted by the loaders (degrees C to K, cm/s to m s-1,
seconds to datenums) are compared to 1e-9 in their units, or to 1 ms for
times; everything else exactly.
"""

import os
import warnings

import netCDF4
import numpy as np
import pytest

from atmosgw.curate.data import availability, problems, states, validate
from atmosgw.curate.obs import (
    load_limb_sounders,
    load_nadir_sounders,
    load_trajectories,
)
from atmosgw.curate.obs._limb_aeolus import profile_groups
from atmosgw.curate.obs._limb_lims import read_lims_file
from atmosgw.curate.obs._traj_dwd_sondes import _parse, read_rows

# ------------------------------------------------------------------- LIMS

HEADER = "%5d%4d%4d%9.4f%9.4f%5d%4d%4d%3d:%2d:%2d%6d%3d%10.4f%6.2f%6.2f%3d%3d"
REST = ("  387100488  -28130831    3370145     0.8588     0.7113   -2.1  348.1"
        "  166.6 -2.803E-05  5  2  2  2  2  2")
CHANNEL = ("   1270   0 2.1562    1      0.000   0      0.000   0      0.000   0"
           "  2.457E-03")


def lims_header(lat, lon, irec, iday, hour, minute, second, idn=2, icloud=100):
    """A header line as the files write it (gt2 and icloud run together when
    icloud has three digits; blanks inside the time)."""
    return (HEADER % (3, 5, 6, lat, lon, 204, irec, iday, hour, minute, second, idn, 0,
                      947.6079, 0.0, 0.0, icloud, 2)) + REST


def write_lims(folder, stamp, year, scans):
    """A LIMS file with 3 levels per scan; scans: list of (header line,
    [3 x 18] level values)."""
    folder.mkdir(parents=True, exist_ok=True)
    lines = ["####  Platform = Nimbus-7", "####  Instrument = LIMS",
             f"####  Year = {year}", "####  Version = V6",
             f"####  Total_Scans = {len(scans)}"]
    for k, (header, levels) in enumerate(scans, start=1):
        lines += ["+" + "#" * 51, f"*###  SCAN {k}", "+" + "#" * 51, header]
        lines += [CHANNEL] * 6
        lines += [" ".join(f"{v:13.6E}" for v in row) for row in levels]
    (folder / f"LIMS_2_{stamp}_V6.ASCII").write_text("\n".join(lines) + "\n")


def lims_levels(temps, heights=(30.0, 25.0, 20.0)):
    rows = np.full((3, 18), 1e24)
    rows[:, 1] = heights            # zt_std
    rows[:, 2] = (10.0, 25.0, 55.0)  # pt_std
    rows[:, 3] = temps              # tt_std
    rows[:, 17] = np.asarray(heights) - 0.3  # geo_heights
    return rows


def test_lims_parsing_and_load(tmp_path):
    # day 312 (1978-11-08): one scan with a glued header, one after midnight
    write_lims(tmp_path / "1978", "19781108", 1978, [
        (lims_header(24.1859, 335.4459, 1, 312, 0, 36, 12),
         lims_levels((226.0, 220.0, 215.0))),
        (lims_header(-55.9227, 269.2260, 2, 312, 2, 4, 5, icloud=90),
         lims_levels((1e24, 221.0, 214.0))),
        (lims_header(10.0, 10.0, 3, 313, 0, 30, 0), lims_levels((227.0, 222.0, 216.0))),
    ])
    scans = read_lims_file(tmp_path / "1978" / "LIMS_2_19781108_V6.ASCII")
    assert scans["n_scans"] == 3 and scans["n_levels"] == 3
    assert scans["alat"][0] == 24.1859 and scans["alon"][1] == 269.2260
    assert scans["icloud"][0] == 100 and scans["second"][1] == 5
    assert np.isnan(scans["tt_std"][1, 0]) and scans["tt_std"][1, 1] == 221.0
    # the 8 Nov file: profiles on 8 Nov from it; the one after midnight is the 9th's
    ds = load_limb_sounders("1978-11-08", "LIMS", file_path=tmp_path, original_z=True,
                            additional_vars=["geo_heights"])
    validate(ds, "profile")
    assert ds.sizes["profile"] == 2
    assert str(ds.time.values[0, 0]) == "1978-11-08T00:36:12.000000"
    assert ds.longitude.values[0, 0] == pytest.approx(335.4459 - 360)
    assert np.allclose(ds.geo_heights.values[0], np.array([30.0, 25.0, 20.0]) - 0.3)
    nxt = load_limb_sounders("1978-11-09", "LIMS", file_path=tmp_path, original_z=True,
                             date_warning=False)
    assert nxt.sizes["profile"] == 1  # (read from the previous day's file)
    assert "LIMS_2_19781108_V6.ASCII" in nxt.attrs["files_read"]
    assert nxt.attrs["files_not_found"] == "1978-11-09"


def test_lims_year_rollover(tmp_path):
    write_lims(tmp_path, "19781231", 1978, [
        (lims_header(0.0, 0.0, 1, 365, 23, 50, 0), lims_levels((220.0, 220.0, 220.0))),
        (lims_header(0.0, 0.0, 2, 1, 0, 10, 0), lims_levels((220.0, 220.0, 220.0))),
    ])
    scans = read_lims_file(tmp_path / "LIMS_2_19781231_V6.ASCII")
    from atmosgw.curate.matlab_time import datenum_to_datetime
    assert str(datenum_to_datetime(scans["time"][1]))[:16] == "1979-01-01T00:10"


def test_lims_before_1979_allowed_only_for_lims(tmp_path):
    with pytest.raises(ValueError, match="satellite era"):
        load_limb_sounders("1978-11-08", "SABER", file_path=tmp_path)


LIMS_FILE = "/data3/LIMS/LIMS_2_19781108_V6.ASCII"


@pytest.mark.skipif(not os.path.isfile(LIMS_FILE), reason="LIMS file absent")
def test_lims_real_matches_users_guide_example():
    # LIMS V6 User's Guide, section 5, "Data Example (day 312)"
    scans = read_lims_file(LIMS_FILE)
    assert (scans["alat"][0], scans["alon"][0]) == (24.1859, 335.4459)
    assert (scans["iorbit"][0], scans["iday"][0], scans["idn"][0]) == (204, 312, 2)
    k = int(np.nanargmin(np.abs(scans["pt_std"][0] - 1.6681)))
    assert scans["zt_std"][0, k] == pytest.approx(44.401, abs=1e-9)
    assert scans["tt_std"][0, k] == pytest.approx(266.908, abs=1e-9)
    assert scans["geo_heights"][0, k] == pytest.approx(44.09225, abs=1e-9)
    # the header format of the synthetic files reproduces a real header
    with open(LIMS_FILE) as handle:
        real = handle.readlines()[11]
    made = lims_header(24.1859, 335.4459, 1, 312, 0, 36, 12, icloud=90)
    assert real.rstrip("\n").replace(" 109 ", "   3 ", 1)[:90] == made[:90]


# ------------------------------------------------------------------ Aeolus


def test_profile_groups_rules():
    # top-down column, a jump back to the top, a repeated bin, a 12 s span
    seconds = np.array([0, 0, 0, 1, 1, 2, 2, 20.0])
    z_bot = np.array([2, 1, 0, 2, 1, 1, 0, 0.0])
    z_top = z_bot + 1
    assert list(profile_groups(seconds, z_bot, z_top)) == [0, 0, 0, 1, 1, 2, 2, 3]


def write_aeolus(folder, duplicate=False):
    """Monthly file 2019_01 with Rayleigh and Mie observations on 2019-01-10
    (duplicate: the first column's three observations twice)."""
    folder.mkdir(parents=True, exist_ok=True)
    t0 = (np.datetime64("2019-01-10T06:00") - np.datetime64("2000-01-01")) / \
        np.timedelta64(1, "s")
    ray = {  # columns top-down; a cloudy (type 1) observation that is not used
        "time": t0 + np.array([0, 0, 0, 0.5, 4, 4, 4]),
        "z_bot": np.array([2000, 1000, 0, 500, 2000, 1000, 0]),
        "hlos": np.array([1000, -250, 32767, 999, 500, 600, 700]),   # cm/s
        "hlos_err": np.array([300, 400, 500, 1, 65535, 200, 200]),
        "obs_type": np.array([2, 2, 2, 1, 2, 2, 2]),
        "valid": np.array([1, 1, 1, 1, 1, 0, 1]),
    }
    if duplicate:  # (copies right after the originals, as in the files)
        ray = {k: np.r_[v[:3], v[:3], v[3:]] for k, v in ray.items()}
    n = ray["time"].size
    with netCDF4.Dataset(folder / "Aeolus_L2B_monthly_2019_01.nc", "w") as f:
        f.Sources = "[('AE_OPER_ALD_U_N_2B_x', '2B16', '03.90')]"
        for ch in ("ray", "mie"):
            f.createDimension(f"{ch}_obs", n)
            for name, dtype in (("time", "f8"), ("lat", "f8"), ("lon", "f8"),
                                ("z_bot", "i4"), ("z_top", "i4"), ("hlos", "i2"),
                                ("hlos_err", "u2"), ("obs_type", "u1"),
                                ("valid", "u1"), ("az", "f8")):
                v = f.createVariable(f"{ch}_{name}", dtype, (f"{ch}_obs",))
                if name in ray:
                    v[:] = ray[name]
                elif name == "z_top":
                    v[:] = ray["z_bot"] + 1000
                elif name == "lat":
                    v[:] = 10 + np.arange(n) * 0.01
                elif name == "lon":
                    v[:] = 20.0
                elif name == "az":
                    v[:] = 260.0
    return ray


def test_aeolus_profiles_units_and_states(tmp_path):
    write_aeolus(tmp_path)
    ds = load_limb_sounders("2019-01-10", "Aeolus", file_path=tmp_path, original_z=True,
                            additional_vars=["z_bot"])
    validate(ds, "profile")
    assert ds.attrs["product"] == "hlos_wind" and ds.attrs["format_version"] == 3
    assert ds.attrs["aeolus_channel"] == "rayleigh"
    assert ds.sizes["profile"] == 2 and ds.sizes["level"] == 3
    assert "temperature" not in ds
    # bottom first; cm/s to m/s; the 32767 fill is missing
    assert np.allclose(ds.altitude.values[0], [0.5, 1.5, 2.5])
    hlos = ds.hlos_wind.values
    assert np.isnan(hlos[0, 0]) and hlos[0, 1] == pytest.approx(-2.5)
    assert hlos[0, 2] == pytest.approx(10.0)
    assert ds.hlos_wind_uncertainty.values[0, 2] == pytest.approx(3.0)
    # the invalid observation (valid 0) is instrument_flag and removed
    flag = ds.hlos_wind_flag.values
    assert states.has_state(flag, "instrument_flag")[1, 1] and np.isnan(hlos[1, 1])
    # uncertainty NaN wherever the wind is
    assert np.all(np.isnan(ds.hlos_wind_uncertainty.values[np.isnan(hlos)]))
    assert ds.hlos_wind.attrs["sign_convention"].startswith("hlos_wind = -u sin(az)")
    assert ds.pressure.attrs["origin"] == "derived"
    assert "L2B baseline 2B16 processor 03.90" in ds.attrs["source_versions"]
    # the cloudy observation is not used by the Rayleigh channel
    assert 999 / 100 not in hlos


def test_aeolus_channel_option(tmp_path):
    write_aeolus(tmp_path)
    ds = load_limb_sounders("2019-01-10", "Aeolus", file_path=tmp_path, original_z=True,
                            aeolus_channel="Mie")
    assert ds.attrs["aeolus_channel"] == "mie"
    assert ds.sizes["profile"] == 1  # (only the type 1 observation)
    with pytest.raises(ValueError, match="aeolus_channel"):
        load_limb_sounders("2019-01-10", "Aeolus", file_path=tmp_path,
                           aeolus_channel="x")


def test_aeolus_previous_month_file(tmp_path):
    # the 2018_12 "monthly" file runs into January: a day in January is found there
    write_aeolus(tmp_path)
    os.rename(tmp_path / "Aeolus_L2B_monthly_2019_01.nc",
              tmp_path / "Aeolus_L2B_monthly_2018_12.nc")
    ds = load_limb_sounders("2019-01-10", "Aeolus", file_path=tmp_path, original_z=True,
                            date_warning=False)
    assert ds.sizes["profile"] == 2
    assert ds.attrs["files_read"] == "Aeolus_L2B_monthly_2018_12.nc"


# ------------------------------------------------------------------- ROMEX


def write_romex(folder):
    folder.mkdir(parents=True, exist_ok=True)
    day = (np.datetime64("2022-10-15") - np.datetime64("1970-01-01")) / \
        np.timedelta64(1, "s")
    with netCDF4.Dataset(folder / "romex_testrmx_100m_v1a.nc", "w") as f:
        f.createDimension("profile", 3)
        f.createDimension("altitude", 4)
        f.createVariable("altitude", "f8", ("altitude",))[:] = [20.0, 20.1, 20.2, 20.3]
        f.createVariable("t", "i8", ("profile",))[:] = [day + 3600, 0, day + 7200]
        f.createVariable("is_bad", "i8", ("profile",))[:] = [0, 0, 1]
        f.createVariable("lat", "f8", ("profile",))[:] = [10.0, 20.0, 30.0]
        f.createVariable("lon", "f8", ("profile",))[:] = [100.0, 110.0, 120.0]
        lat = np.tile([[10.0], [20.0], [30.0]], (1, 4))
        lat[0, 3] = np.nan
        f.createVariable("lat_full", "f8", ("profile", "altitude"))[:] = lat
        f.createVariable("lon_full", "f8", ("profile", "altitude"))[:] = \
            np.tile([[100.0], [110.0], [120.0]], (1, 4))
        f.createVariable("temp", "f8", ("profile", "altitude"))[:] = -60.0
        f.createVariable("press", "f8", ("profile", "altitude"))[:] = 50.0


def test_romex(tmp_path):
    write_romex(tmp_path)
    ds = load_limb_sounders("2022-10-15", "ROMEX", file_path=tmp_path, original_z=True)
    validate(ds, "profile")
    assert ds.sizes["profile"] == 2  # (the profile with t = 0 is not returned)
    assert ds.temperature.values[0, 0] == pytest.approx(213.15, abs=1e-9)
    assert ds.latitude.values[0, 3] == 10.0  # (lat_full NaN: the profile's lat)
    flag = ds.temperature_flag.values
    assert np.all(states.has_state(flag[1], "instrument_flag"))  # is_bad
    assert np.all(np.isnan(ds.temperature.values[1]))
    assert "profiles with t = 0 are not returned (1" in ds.time.attrs["comment"]
    assert ds.attrs["source_versions"].startswith("ROMEX v1a")
    kept = load_limb_sounders("2022-10-15", "ROMEX", file_path=tmp_path,
                              original_z=True, instrument_qc=False)
    assert np.all(np.isfinite(kept.temperature.values[1]))


# --------------------------------------------------- sampled data (SWANS)


def test_sampler_netcdf_and_sampled_from(tmp_path):
    n_lev, n_prof = 4, 3
    with netCDF4.Dataset(tmp_path / "HIRDLS-test_model-2020d020.nc", "w") as f:
        f.createDimension("dim_3", n_lev)
        f.createDimension("dim_4", n_prof)
        def var(name, value):
            v = f.createVariable(name, "f8", ("dim_3", "dim_4"), fill_value=-9999.0)
            v[:] = value
        alt = np.tile(np.array([20.0, 25, 30, 35])[:, None], (1, n_prof))
        var("Alt", alt)
        var("Lat", np.tile(np.array([1.0, 2, 3]), (n_lev, 1)))
        var("Lon", np.tile(np.array([4.0, 5, 6]), (n_lev, 1)))
        day = 737810.0  # 2020-01-20
        var("Time", np.tile(day + np.array([0.1, 0.2, 0.3]), (n_lev, 1)))
        var("Prs", 1013.25 * np.exp(-alt / 7))
        temp = 200 + alt + np.arange(n_prof)[None, :]
        temp[0, 2] = -9999.0
        var("Temp", temp)
        var("ObsProf", np.tile(np.array([1.0, 2, 3]), (n_lev, 1)))
    ds = load_limb_sounders("2020-01-20", "Sampler", sample_path=tmp_path,
                            original_z=True, sampled_from="test model",
                            file_source=True)
    validate(ds, "profile")
    assert ds.sizes["profile"] == 3
    assert ds.temperature.values[1, 2] == 231.0  # (profile 2, 30 km)
    assert np.isnan(ds.temperature.values[2, 0])  # (-9999 fill)
    assert ds.latitude.values[2, 0] == 3.0
    assert list(ds.source_profile.values[:, 0]) == [1, 2, 3]
    assert ds.attrs["data"] == "sampled" and ds.attrs["sampled_from"] == "test model"
    assert ds.temperature.attrs["origin"] == "derived"


def test_nadir_sampled_from(tmp_path):
    from test_nadir_known import D, seconds, write_granule
    write_granule(tmp_path, 53, seconds(53))
    plain = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path)
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                             sampled_from="a model")
    validate(ds, "swath")
    assert plain.temperature.attrs["origin"] == "retrieved"
    assert plain.attrs["data"] == "raw"
    assert ds.temperature.attrs["origin"] == "derived" and ds.attrs["data"] == "sampled"
    assert ds.attrs["sampled_from"] == "a model"
    finite = np.isfinite(ds.temperature.values)
    assert np.all(states.has_state(ds.temperature_flag.values, "derived") == finite)
    np.testing.assert_array_equal(ds.temperature.values, plain.temperature.values)
    with pytest.raises(TypeError):
        load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                            sampled_from=3)


# ------------------------------------------------------- DWD radiosondes

PRODUCT_HEADER = ("STATIONS_ID;BEZUGSDATUM_SYNOP;MESSZEITPUNKT;QN_1;AE_GB_POS;"
                  "AE_GL_POS;AE_GPM;AE_P;AE_TT;AE_RF;AE_TD;AE_FF;AE_DD;eor")


def write_sondes(root, n_days=31, samples=60):
    """Station 01234 (Teststadt, 50 N 10 E), January 2020: ascents at 05 and
    11 (reference hours), launched 15 min before; the 2020-01-03T05 ascent
    has no positions; one sample has latitude 0.00 (a fill value)."""
    folder = root / "2020"
    folder.mkdir(parents=True)
    lines, meta = [PRODUCT_HEADER], [
        "Stations_id;Bezugsdatum;AE_Start;AE_Ende;Gipfel_hpa;Gipfel_gpm;P0-Boden;"
        "TT-Boden;RF-Boden;Geraete_Name;Geraete_Art;eor"]
    for day in range(1, n_days + 1):
        for hour in (5, 11):
            ref = f"202001{day:02d}{hour:02d}"
            meta.append(f"  1234;{ref};202001{day:02d}{hour - 1:02d}45;"
                        f"202001{day:02d}{hour:02d}45;  10; 30000; 1000.0;  1.0;"
                        "  80.0;Radiosonde TEST;MSGRT;eor;")
            for k in range(samples):
                no_pos = day == 3 and hour == 5
                lat = ("-999" if no_pos else "0.00" if k == 7
                       else f"{50 + k * 0.001:.2f}")
                lon = "-999" if no_pos else f"{10 + k * 0.01:.2f}"
                temp = "-223.2" if (day == 2 and k == 5) else f"{10 - k:.1f}"
                lines.append(f"  1234;{ref};{k * 2:6d}; -999; {lat} ; {lon} ; "
                             f"{100 + 10 * k}; {1000 - k:.2f}; {temp} ; 80.0 ; "
                             f"0.0 ; 10.0 ; 270;eor")
    (folder / "produkt_sec_aero_20200101_20201231_01234.txt").write_text(
        "\n".join(lines) + "\n")
    (folder / "Metadaten_Sekunde_Aero_01234.txt").write_text("\n".join(meta) + "\n",
                                                             encoding="latin-1")
    (folder / "Metadaten_Geographie_01234.txt").write_text(
        "Stations_id;Stationshoehe;Hoehe_Messfeld_Aerologie;Geogr.Breite;"
        "Geogr.Laenge;von_datum;bis_datum;Stationsname\n"
        "  1234;  100.00;        ; 50.0000; 10.0000;19900101;        ;Teststadt\n",
        encoding="latin-1")
    return folder / "produkt_sec_aero_20200101_20201231_01234.txt"


def test_sonde_bisection_equals_full_read(tmp_path):
    path = write_sondes(tmp_path)
    assert os.path.getsize(path) > 4 * 65536  # (so the bisection is used)
    with open(path, "rb") as handle:
        handle.readline()
        every = _parse(handle.readlines())
    for lo, hi in ((0, 2020010223), (2020011500, 2020011605), (2020013100, 9999999999),
                   (2020020100, 2020020223)):
        rows = read_rows(path, lo, hi)
        expected = every[(every[:, 1] >= lo) & (every[:, 1] <= hi)]
        np.testing.assert_array_equal(rows, expected)


def test_sondes_load(tmp_path):
    write_sondes(tmp_path)
    ds = load_trajectories("2020-01-02", "DWDsondes", file_path=tmp_path,
                           additional_vars=["AE_RF"], file_source=True)
    validate(ds, "trajectory")
    assert ds.sizes["flight"] == 2 and ds.sizes["sample"] == 60
    assert list(ds.flight_id.values) == ["01234-2020-01-02T04:45",
                                         "01234-2020-01-02T10:45"]
    assert list(ds.platform_name.values) == ["Teststadt"] * 2
    assert str(ds.start_time.values[0]) == "2020-01-02T04:45:00.000000"
    assert str(ds.time.values[0, 3]) == "2020-01-02T04:45:06.000000"  # (+ 2 s x 3)
    assert ds.altitude.values[0, 2] == pytest.approx(0.12)
    assert ds.temperature.values[0, 0] == pytest.approx(283.15, abs=1e-9)
    # wind from the west (270 degrees) at 10 m/s: eastward +10
    assert ds.eastward_wind.values[0, 0] == pytest.approx(10.0)
    assert abs(ds.northward_wind.values[0, 0]) < 1e-9
    # latitude exactly 0.00 is a fill value
    assert np.isnan(ds.latitude.values[0, 7]) and ds.latitude.values[0, 8] == 50.01
    # 50 K: outside the physical range, removed from temperature only
    flag = ds.temperature_flag.values
    assert states.has_state(flag, "outside_physical_range")[0, 5]
    assert np.isnan(ds.temperature.values[0, 5])
    assert np.isfinite(ds.pressure.values[0, 5])
    assert ds.temperature.attrs["origin"] == "measured"
    assert ds.AE_RF.values[0, 0] == 80.0 and ds.original_files.size == 1


def test_sondes_station_position_and_filters(tmp_path):
    write_sondes(tmp_path)
    ds = load_trajectories("2020-01-03", "DWDsondes", file_path=tmp_path,
                           stations="teststadt")
    assert ds.sizes["flight"] == 2
    assert np.all(ds.latitude.values[0] == 50.0)
    assert np.all(ds.longitude.values[0] == 10.0)
    assert "(1 of 2 flights)" in ds.latitude.attrs["comment"]
    assert load_trajectories("2020-01-03", "DWDsondes", file_path=tmp_path,
                             stations=1234).sizes["flight"] == 2
    with pytest.warns(UserWarning):
        none = load_trajectories("2020-01-03", "DWDsondes", file_path=tmp_path,
                                 stations=["99999"])
    assert none.sizes["flight"] == 0
    assert "station 99999" in none.attrs["files_not_found"]
    validate(none, "trajectory")


def test_sonde_launch_before_midnight_belongs_to_that_day(tmp_path):
    # an ascent for 00 UTC launched at 22:45 the day before
    write_sondes(tmp_path, n_days=2)
    folder = tmp_path / "2020"
    text = (folder / "produkt_sec_aero_20200101_20201231_01234.txt").read_text()
    text += ("  1234;2020010300;     0; -999; 50.00 ; 10.00 ; 100; 1000.00; 5.0 ;"
             " 80.0 ; 0.0 ; 1.0 ; 90;eor\n")
    (folder / "produkt_sec_aero_20200101_20201231_01234.txt").write_text(text)
    meta = folder / "Metadaten_Sekunde_Aero_01234.txt"
    with open(meta, "a", encoding="latin-1") as f:
        f.write("  1234;2020010300;202001022245;202001030045;  10; 30000; 1000.0;  1.0;"
                "  80.0;Radiosonde TEST;MSGRT;eor;\n")
    ds = load_trajectories("2020-01-02", "DWDsondes", file_path=tmp_path)
    assert ds.sizes["flight"] == 3
    assert str(ds.start_time.values[-1]) == "2020-01-02T22:45:00.000000"
    assert str(ds.reference_time.values[-1]) == "2020-01-03T00:00:00.000000"


# ----------------------------------------------------------------- formats


def test_profile_format_versions_and_products(tmp_path):
    write_romex(tmp_path)
    ds = load_limb_sounders("2022-10-15", "ROMEX", file_path=tmp_path, original_z=True)
    # version 3 tracks the winds: listed, as not measured
    assert availability(ds, "eastward_wind") == "not_measured"
    old = ds.copy()
    old.attrs["format_version"] = 1
    validate(old, "profile")  # (version 1 Datasets are still accepted)
    # an older Dataset without the version 3 entries is still valid as itself
    listed = states.parse_availability(ds.attrs["availability"])
    for name in ("eastward_wind", "northward_wind", "upward_wind"):
        listed.pop(name)
    older = ds.copy()
    older.attrs.update(format_version=2,
                       availability=states.format_availability(listed))
    validate(older, "profile")
    older.attrs["format_version"] = 3
    assert any("eastward_wind" in p for p in problems(older, "profile"))
    bad = ds.copy()
    bad.attrs["format_version"] = 4
    assert any("format_version" in p for p in problems(bad, "profile"))
    wind = ds.copy()
    wind.attrs["product"] = "hlos_wind"
    found = problems(wind, "profile")
    assert any("hlos_wind" in p for p in found)


def test_trajectory_needs_a_measured_variable(tmp_path):
    write_sondes(tmp_path, n_days=2)
    ds = load_trajectories("2020-01-02", "DWDsondes", file_path=tmp_path)
    stripped = ds.drop_vars([v for v in ds.data_vars if v.split("_flag")[0] in
                             ("temperature", "pressure", "eastward_wind",
                              "northward_wind")])
    assert any("needs at least one of" in p for p in problems(stripped, "trajectory"))


# ------------------------------------------- fixes after the review (2026-11-27)


def test_lims_profiles_from_the_next_days_file(tmp_path):
    # the 27 March file starts before 27 March: its first scan belongs to the 26th
    write_lims(tmp_path, "19790326", 1979, [
        (lims_header(0.0, 0.0, 1, 85, 1, 0, 0), lims_levels((220.0, 220.0, 220.0)))])
    write_lims(tmp_path, "19790327", 1979, [
        (lims_header(0.0, 0.0, 1, 85, 23, 0, 0), lims_levels((221.0, 221.0, 221.0))),
        (lims_header(0.0, 0.0, 2, 86, 1, 0, 0), lims_levels((222.0, 222.0, 222.0)))])
    ds = load_limb_sounders("1979-03-26", "LIMS", file_path=tmp_path, original_z=True,
                            date_warning=False)
    assert ds.sizes["profile"] == 2
    assert list(ds.temperature.values[:, 0]) == [220.0, 221.0]


def test_aeolus_duplicates_read_once(tmp_path):
    write_aeolus(tmp_path / "plain")
    write_aeolus(tmp_path / "copies", duplicate=True)
    plain, once = (load_limb_sounders("2019-01-10", "Aeolus", file_path=tmp_path / k,
                                      original_z=True) for k in ("plain", "copies"))
    assert once.sizes == plain.sizes
    np.testing.assert_array_equal(once.hlos_wind.values, plain.hlos_wind.values)


def test_sondes_copies_of_a_file_read_once(tmp_path):
    path = write_sondes(tmp_path, n_days=3)
    copy = tmp_path / "1995"
    copy.mkdir()
    for name in os.listdir(path.parent):
        (copy / name).write_bytes((path.parent / name).read_bytes()[:-2000])
    with pytest.warns(UserWarning, match="more than one"):
        ds = load_trajectories("2020-01-02", "DWDsondes", file_path=tmp_path,
                               file_source=True)
    assert ds.sizes["flight"] == 2 and ds.original_files.size == 1
    assert "used 2020" in ds.attrs["multiple_matches"]
    validate(ds, "trajectory")
    repeated = ds.copy()
    repeated["flight_id"] = ("flight", np.array(["a", "a"]))
    assert any("repeated" in p for p in problems(repeated, "trajectory"))


def test_hlos_product_needs_version_2(tmp_path):
    write_aeolus(tmp_path)
    ds = load_limb_sounders("2019-01-10", "Aeolus", file_path=tmp_path, original_z=True)
    ds.attrs["format_version"] = 1
    assert any("version 2" in p for p in problems(ds, "profile"))
    with pytest.warns(UserWarning, match="below 18 km"):
        load_limb_sounders("2019-01-10", "Aeolus", file_path=tmp_path)


# ------------------------------------------------------------ IAGOS, FAAM


def iagos_file(path, new=True, n=50):
    """An IAGOS flight departing 2019-03-18T16:40:00 (new: L2 3.1 style with
    validity flags and a package temperature; old: MOZAIC style)."""
    kind = "NETCDF4" if new else "NETCDF3_CLASSIC"
    with netCDF4.Dataset(path, "w", format=kind) as f:
        f.createDimension("UTC_time", n)
        t = f.createVariable("UTC_time", "f8", ("UTC_time",),
                             fill_value=0.0 if new else -9999.0)
        t.units = "seconds since 2019-03-18 00:00:00" + (" UTC" if new else "")
        t[:] = 60000 + 4.0 * np.arange(n)
        def var(name, value, flag=None):
            v = f.createVariable(name, "f8", ("UTC_time",), fill_value=-9999.0)
            v[:] = value
            if flag is not None and new:
                fv = f.createVariable(f"{name}_validity_flag", "i4", ("UTC_time",),
                                      fill_value=-9999)
                fv[:] = flag
        flags = np.zeros(n, dtype=int)
        flags[3], flags[4], flags[5] = 3, 2, 7  # erroneous, limited, missing
        var("lat", np.linspace(50, 52, n))
        var("lon", np.linspace(8, 10, n))
        var("baro_alt_AC", np.full(n, 10000.0))
        var("gps_alt_AC", np.full(n, 10100.0) if new else np.full(n, -9999.0))
        pres = np.full(n, 26500.0)
        pres[7] = -999900.0
        var("air_press_AC", pres, flags)
        var("air_temp_AC", np.full(n, 220.0), flags)
        var("air_temp_P1" if new else "air_temp_PM", np.full(n, 221.0), flags)
        var("zon_wind_AC", np.full(n, 30.0), flags)
        var("mer_wind_AC", np.full(n, -5.0), flags)
        f.flight_name = "2019031816400001"
        f.departure_UTC_time = "2019-03-18T16:40:00Z"
        f.departure_airport = "FRA, Frankfurt, Germany"
        f.arrival_airport = "CPT, Cape Town, South Africa"
        if new:
            f.platform = "Aircraft, IAGOS-01, MSN-304, D-AIGT, Lufthansa"
            f.file_version = "3.1.2"
        else:
            f.title = "Instrument-75, MOZAIC-05, A340, OE-LAG, Austrian"
            f.source = "IAGOS-MOZAIC"
            f.file_revision = "1.0 20191021"


def test_iagos_from_zip_and_folder(tmp_path):
    import zipfile
    work = tmp_path / "work"
    work.mkdir()
    iagos_file(work / "IAGOS_timeseries_2019031816400001_L2_3.1.2.nc4")
    with zipfile.ZipFile(tmp_path / "iagos.zip", "w") as archive:
        archive.write(work / "IAGOS_timeseries_2019031816400001_L2_3.1.2.nc4",
                      "TimeSeries/iagos_2019/IAGOS_timeseries_2019031816400001_L2_3.1.2.nc4")
    (tmp_path / "zip").mkdir()
    (tmp_path / "iagos.zip").rename(tmp_path / "zip" / "iagos.zip")
    from_zip = load_trajectories("2019-03-18", "IAGOS", file_path=tmp_path / "zip",
                                 additional_vars=["air_temp_AC"])
    from_folder = load_trajectories("2019-03-18", "IAGOS", file_path=work)
    with warnings.catch_warnings():  # (zip and an unzipped copy: no warning)
        warnings.simplefilter("error")
        both = load_trajectories("2019-03-18", "IAGOS", file_path=tmp_path)
    assert both.sizes["flight"] == 1 and both.attrs["multiple_matches"] == ""
    validate(from_zip, "trajectory")
    np.testing.assert_array_equal(from_zip.temperature.values,
                                  from_folder.temperature.values)
    ds = from_zip
    assert ds.sizes["flight"] == 1 and ds.flight_id.values[0] == "2019031816400001"
    assert ds.platform_id.values[0] == "D-AIGT"
    assert ds.temperature_source.values[0] == "air_temp_P1"
    assert ds.altitude_source.values[0] == "gps_alt_AC"
    assert ds.altitude.values[0, 0] == pytest.approx(10.1)
    assert ds.pressure.values[0, 0] == pytest.approx(265.0)
    assert np.isnan(ds.pressure.values[0, 7])  # (-999900 Pa)
    t = ds.temperature.values[0]
    assert t[0] == 221.0 and t[4] == 221.0  # (limited: kept)
    assert np.isnan(t[3]) and np.isnan(t[5])  # (erroneous, missing: removed)
    assert states.has_state(ds.temperature_flag.values[0], "instrument_flag")[3]
    assert str(ds.time.values[0, 1]) == "2019-03-18T16:40:04.000000"
    assert ds.air_temp_AC.values[0, 0] == 220.0
    assert ds.departure_airport.values[0].startswith("FRA")
    none = load_trajectories("2019-03-20", "IAGOS", file_path=tmp_path / "zip",
                             date_warning=False)
    assert none.sizes["flight"] == 0


def test_iagos_copies_zip_backup_and_versions(tmp_path):
    import zipfile
    work = tmp_path / "work"
    work.mkdir()
    a = work / "IAGOS_timeseries_2019031816400001_L2_3.1.2.nc4"
    b = work / "IAGOS_timeseries_2019031816400002_L2_3.1.2.nc4"
    iagos_file(a)
    iagos_file(b)
    with netCDF4.Dataset(b, "a") as f:
        f.flight_name = "2019031816400002"
    data = tmp_path / "IAGOS"
    data.mkdir()
    with zipfile.ZipFile(data / "iagos.zip", "w") as archive:
        for path in (a, b):
            archive.write(path, f"TimeSeries/iagos_2019/{path.name}")
    # flight 01 unzipped beside the archive as a newer version with other
    # values; flight 02 only in the archive
    unzipped = data / "TimeSeries" / "iagos_2019"
    unzipped.mkdir(parents=True)
    newer = unzipped / "IAGOS_timeseries_2019031816400001_L2_3.2.0.nc4"
    iagos_file(newer)
    with netCDF4.Dataset(newer, "a") as f:
        f["air_temp_P1"][0] = 230.0
        f.file_version = "3.2.0"
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ds = load_trajectories("2019-03-18", "IAGOS", file_path=data)
    assert ds.sizes["flight"] == 2
    assert sorted(ds.flight_id.values) == ["2019031816400001", "2019031816400002"]
    one = int(np.flatnonzero(ds.flight_id.values == "2019031816400001")[0])
    assert ds.temperature.values[one, 0] == 230.0  # (the unzipped copy)
    assert ds.attrs["multiple_matches"] == ""
    # two unzipped versions of a flight: the newest, recorded
    older = unzipped / "IAGOS_timeseries_2019031816400001_L2_3.1.2.nc4"
    iagos_file(older)
    with pytest.warns(UserWarning, match="more than one"):
        ds = load_trajectories("2019-03-18", "IAGOS", file_path=data)
    one = int(np.flatnonzero(ds.flight_id.values == "2019031816400001")[0])
    assert ds.temperature.values[one, 0] == 230.0
    assert "3.2.0" in ds.attrs["multiple_matches"]


def test_iagos_old_files(tmp_path):
    iagos_file(tmp_path / "IAGOS_timeseries_2019031816400001.nc", new=False)
    ds = load_trajectories("2019-03-18", "IAGOS", file_path=tmp_path)
    validate(ds, "trajectory")
    assert ds.temperature_source.values[0] == "air_temp_PM"
    assert ds.altitude_source.values[0] == "baro_alt_AC"  # (no GPS altitude)
    assert ds.altitude.values[0, 0] == pytest.approx(10.0)
    assert ds.platform_id.values[0] == "OE-LAG"
    assert np.all(np.isfinite(ds.temperature.values[0]))  # (no flags: no screening)


def faam_file(path, n=4):
    with netCDF4.Dataset(path, "w") as f:
        f.createDimension("Time", n)
        f.createDimension("sps32", 32)
        t = f.createVariable("Time", "i4", ("Time",), fill_value=-9999)
        t.units = "seconds since 2025-07-09 00:00:00 +0000"
        t[:] = 36000 + np.arange(n)
        def var(name, value, flagged=None):
            v = f.createVariable(name, "f4", ("Time", "sps32"), fill_value=-9999.0)
            v[:] = value
            fl = f.createVariable(f"{name}_FLAG", "i1", ("Time", "sps32"))
            fl[:] = 0 if flagged is None else flagged
        ramp = np.tile(np.arange(32.0), (n, 1))
        var("LAT_GIN", 50 + 0 * ramp)
        var("LON_GIN", -1 + 0 * ramp)
        var("ALT_GIN", 5000 + 0 * ramp)
        var("PALT_RVS", 4900 + 0 * ramp)
        var("PS_RVSM", 540 + 0 * ramp)
        flagged = np.zeros((n, 32), dtype=int)
        flagged[1, :16] = 4  # (in cloud: half of second 1)
        flagged[2, :] = 1  # (all of second 2)
        var("TAT_DI_R", 260 + ramp, flagged)
        var("U_C", 10 + 0 * ramp)
        var("V_C", 0 * ramp)
        var("W_C", 0.1 + 0 * ramp)
        f.flight_number = "c999"
        f.platform = "FAAM BAe-146-301 Atmospheric Research Aircraft"
        f.time_coverage_start = "2025-07-09T10:00:00Z"
        f.revision_number = "0"
        f.processing_software_version = "25.0.0"


def test_faam_rates_and_flags(tmp_path):
    faam_file(tmp_path / "core_faam_20250709_v005_r0_c999.nc")
    ds = load_trajectories("2025-07-09", "FAAM", file_path=tmp_path,
                           additional_vars=["PALT_RVS"])
    validate(ds, "trajectory")
    assert ds.attrs["sample_rate"] == "1hz" and ds.sizes["sample"] == 4
    t = ds.temperature.values[0]
    assert t[0] == pytest.approx(260 + 15.5)  # (mean of 0..31)
    assert t[1] == pytest.approx(260 + 23.5)  # (mean of the unflagged 16..31)
    assert np.isnan(t[2])
    assert states.has_state(ds.temperature_flag.values[0], "instrument_flag")[2]
    assert ds.altitude.values[0, 0] == pytest.approx(5.0)
    assert ds.upward_wind.values[0, 0] == pytest.approx(0.1, abs=1e-6)
    assert str(ds.time.values[0, 0])[:23] == "2025-07-09T10:00:00.484"
    full = load_trajectories("2025-07-09", "FAAM", file_path=tmp_path,
                             sample_rate="full")
    validate(full, "trajectory")
    assert full.sizes["sample"] == 128
    assert full.temperature.values[0, 5] == 265.0
    assert np.isnan(full.temperature.values[0, 32 + 3])  # (flagged sample)
    assert np.isfinite(full.temperature.values[0, 32 + 20])
    # (1/32 s; a datenum resolves about 10 us today, hence the tolerance)
    step = np.diff(full.time.values[0, :3]).astype("timedelta64[us]").astype(int)
    assert np.all(np.abs(step - 31250) <= 20)
    with pytest.raises(TypeError, match="sample_rate"):
        load_trajectories("2019-03-18", "IAGOS", file_path=tmp_path, sample_rate="full")
    with pytest.raises(TypeError, match="stations"):
        load_trajectories("2025-07-09", "FAAM", file_path=tmp_path, stations=1)


def test_iagos_review_fixes(tmp_path):
    import zipfile
    # a flight whose package temperature is all flagged erroneous: the
    # aircraft's is used
    work = tmp_path / "work"
    work.mkdir()
    path = work / "IAGOS_timeseries_2019031816400001_L2_3.1.2.nc4"
    iagos_file(path)
    with netCDF4.Dataset(path, "a") as f:
        f["air_temp_P1_validity_flag"][:] = 3
    ds = load_trajectories("2019-03-18", "IAGOS", file_path=work)
    assert ds.temperature_source.values[0] == "air_temp_AC"
    assert ds.temperature.values[0, 0] == 220.0
    # copies of a flight on another day are not reported for this day
    other = work / "IAGOS_timeseries_2019010100000001_L2_3.1.2.nc4"
    iagos_file(other)
    with zipfile.ZipFile(tmp_path / "iagos.zip", "w") as archive:
        for folder in ("iagos_2018", "iagos_2019"):
            archive.write(other, f"TimeSeries/{folder}/{other.name}")
        archive.write(path, f"TimeSeries/iagos_2019/{path.name}")
    (tmp_path / "zip").mkdir()
    (tmp_path / "iagos.zip").rename(tmp_path / "zip" / "iagos.zip")
    ds = load_trajectories("2019-03-18", "IAGOS", file_path=tmp_path / "zip")
    assert ds.attrs["multiple_matches"] == ""
    with pytest.raises(TypeError, match="sample_rate"):
        load_trajectories("2019-03-18", "IAGOS", file_path=work, sample_rate="FULL")
    assert load_trajectories("2019-03-18", "IAGOS", file_path=work,
                             sample_rate="1HZ").sizes["flight"] == 1


def test_faam_review_fixes(tmp_path):
    faam_file(tmp_path / "core_faam_20250709_v005_r0_c999.nc")
    with netCDF4.Dataset(tmp_path / "core_faam_20250709_v005_r0_c999.nc", "a") as f:
        f.createDimension("sps20", 20)
        v = f.createVariable("X20", "f4", ("Time", "sps20"), fill_value=-9999.0)
        v[:] = 1.0
        f["LAT_GIN_FLAG"][0, :] = 1  # (no position for second 0)
    with pytest.raises(ValueError, match="20 Hz"):
        load_trajectories("2025-07-09", "FAAM", file_path=tmp_path, sample_rate="full",
                          additional_vars=["X20"])
    ds = load_trajectories("2025-07-09", "FAAM", file_path=tmp_path,
                           additional_vars=["X20"])
    assert ds.X20.values[0, 0] == 1.0  # (1 Hz: any rate)
    assert np.isnan(ds.latitude.values[0, 0]) and ds.latitude.values[0, 1] == 50.0
    full = load_trajectories("2025-07-09", "FAAM", file_path=tmp_path,
                             sample_rate="full")
    assert np.all(np.isnan(full.latitude.values[0, :32]))
    # a second revision of the flight: the highest is read
    faam_file(tmp_path / "core_faam_20250709_v005_r1_c999.nc")
    with pytest.warns(UserWarning, match="more than one"):
        ds = load_trajectories("2025-07-09", "FAAM", file_path=tmp_path,
                               file_source=True)
    assert ds.sizes["flight"] == 1
    assert ds.original_files.values[0] == "core_faam_20250709_v005_r1_c999.nc"


# ------------------------------- versions, availability, records (2026-11-27b)


def _iagos_copy(path, temp, version=None, modified=None):
    """An IAGOS flight file with temperature temp (and its stated dates)."""
    iagos_file(path)
    with netCDF4.Dataset(path, "a") as f:
        f["air_temp_P1"][:] = temp
        if version is not None:
            f.file_version = version
        if modified is not None:
            f.date_modified = modified


def test_iagos_newest_version_then_date_then_unzipped(tmp_path):
    import zipfile
    work = tmp_path / "work"
    work.mkdir()
    data = tmp_path / "IAGOS"
    unzipped = data / "TimeSeries" / "iagos_2019"
    unzipped.mkdir(parents=True)
    stem = "IAGOS_timeseries_2019031816400001_L2_"
    member = f"TimeSeries/iagos_2019/{stem}3.2.0.nc4"
    # (1) a newer version in the zip beats an older unzipped copy (recorded)
    _iagos_copy(work / f"{stem}3.2.0.nc4", 232.0, "3.2.0")
    with zipfile.ZipFile(data / "iagos.zip", "w") as archive:
        archive.write(work / f"{stem}3.2.0.nc4", member)
    _iagos_copy(unzipped / f"{stem}3.1.2.nc4", 231.0)
    with pytest.warns(UserWarning, match="more than one"):
        ds = load_trajectories("2019-03-18", "IAGOS", file_path=data)
    assert ds.temperature.values[0, 0] == 232.0
    assert "3.2.0" in ds.attrs["multiple_matches"]
    # (2) the same version: the later date stated in the file wins, even zipped
    (unzipped / f"{stem}3.1.2.nc4").unlink()
    _iagos_copy(unzipped / f"{stem}3.2.0.nc4", 233.0, "3.2.0", "2023-03-03T14:37:51Z")
    with zipfile.ZipFile(data / "iagos.zip", "w") as archive:
        _iagos_copy(work / f"{stem}3.2.0.nc4", 234.0, "3.2.0", "2024-01-01T00:00:00Z")
        archive.write(work / f"{stem}3.2.0.nc4", member)
    with pytest.warns(UserWarning, match="more than one"):
        ds = load_trajectories("2019-03-18", "IAGOS", file_path=data)
    assert ds.temperature.values[0, 0] == 234.0
    # (3) the same version and date: unzipped first, silently (the archive's
    # backup)
    with zipfile.ZipFile(data / "iagos.zip", "w") as archive:
        _iagos_copy(work / f"{stem}3.2.0.nc4", 235.0, "3.2.0", "2023-03-03T14:37:51Z")
        archive.write(work / f"{stem}3.2.0.nc4", member)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ds = load_trajectories("2019-03-18", "IAGOS", file_path=data)
    assert ds.temperature.values[0, 0] == 233.0 and ds.attrs["multiple_matches"] == ""


def test_newest_version_in_file_names(tmp_path):
    from atmosgw.curate.obs._common import choose_file
    from atmosgw.curate.obs._limb_tools import VERSION_IN_NAME, newest_version
    cases = {
        "HIRDLS": ["HIRDLS-Aura_L2_v06-00-00-c01_2007d001.he5",
                   "HIRDLS-Aura_L2_v07-00-20-c01_2007d001.he5"],
        "SABER": ["SABER_Temp_O3_H2O_January2007_v2.0.nc",
                  "SABER_Temp_O3_H2O_January2007_v2.07.nc"],
        "MIPAS": ["mipas_envisat_20070102_nom_oxf_v1.40.nc",
                  "mipas_envisat_20070102_nom_oxf_v1.5.nc"],
        "SOFIE": ["SOFIE_Level2_2007135_01.2.nc", "SOFIE_Level2_2007135_01.3.nc"],
    }
    for instrument, names in cases.items():
        folder = tmp_path / instrument
        (folder / "a").mkdir(parents=True)
        (folder / "b").mkdir()
        # the newer version is second in search order, so MATLAB would take
        # the older one
        paths = [str(folder / "a" / names[0]), str(folder / "b" / names[1])]
        for path in paths:
            open(path, "w").close()
        info = {"path": str(folder), "record": {"multiple_matches": []}}
        chosen = choose_file(paths, info, "day",
                             newest_version(VERSION_IN_NAME[instrument]))
        newest = paths[0] if instrument == "MIPAS" else paths[1]  # (1.40 > 1.5)
        assert chosen == newest, instrument
        assert len(info["record"]["multiple_matches"]) == 1
        # equal (or no) versions: the first, as MATLAB
        assert choose_file([paths[0], paths[0]], info, "day",
                           newest_version(VERSION_IN_NAME[instrument])) == paths[0]


def test_measured_quantities_always_listed(tmp_path):
    from test_nadir_known import D, seconds, write_granule
    write_granule(tmp_path / "airs", 53, seconds(53))
    swath = load_nadir_sounders(None, "AIRS", granules=[(D, 53)],
                                file_path=tmp_path / "airs")
    assert swath.attrs["format_version"] == 2
    assert availability(swath, "pressure") == "not_in_files"
    for name in ("eastward_wind", "northward_wind", "upward_wind"):
        assert availability(swath, name) == "not_measured"
    write_romex(tmp_path / "romex")
    limb = load_limb_sounders("2022-10-15", "ROMEX", file_path=tmp_path / "romex",
                              original_z=True)
    assert availability(limb, "upward_wind") == "not_measured"
    write_aeolus(tmp_path / "aeolus")
    wind = load_limb_sounders("2019-01-10", "Aeolus", file_path=tmp_path / "aeolus",
                              original_z=True)
    validate(wind, "profile")
    assert wind.sizes["profile"] > 0
    assert "instrument_flag: observations with ray_valid 0" in \
        wind.hlos_wind.attrs["comment"]
    assert availability(wind, "temperature") == "not_measured"
    write_sondes(tmp_path / "sondes", n_days=2)
    sondes = load_trajectories("2020-01-02", "DWDsondes", file_path=tmp_path / "sondes")
    assert availability(sondes, "upward_wind") == "not_measured"
    # an older swath (version 1, no such entries) is still valid
    listed = states.parse_availability(swath.attrs["availability"])
    for name in ("pressure", "eastward_wind", "northward_wind", "upward_wind"):
        listed.pop(name)
    old = swath.copy()
    old.attrs.update(format_version=1, availability=states.format_availability(listed))
    validate(old, "swath")


def test_screening_rules_in_comments_and_empty_files_recorded(tmp_path):
    write_romex(tmp_path / "romex")
    ds = load_limb_sounders("2022-10-15", "ROMEX", file_path=tmp_path / "romex",
                            original_z=True)
    assert "instrument_flag: profiles the file marks is_bad" in \
        ds.temperature.attrs["comment"]
    assert "is_bad" in ds.pressure.attrs["comment"]
    (tmp_path / "iagos").mkdir()
    iagos_file(tmp_path / "iagos" / "IAGOS_timeseries_2019031816400001_L2_3.1.2.nc4")
    flights = load_trajectories("2019-03-18", "IAGOS", file_path=tmp_path / "iagos")
    for name in ("temperature", "pressure", "eastward_wind", "northward_wind"):
        assert "instrument_flag: validity flag 3" in flights[name].attrs["comment"]
    # a ROMEX file whose temperatures for the day are all missing
    with netCDF4.Dataset(tmp_path / "romex" / "romex_testrmx_100m_v1a.nc", "a") as f:
        f["temp"][:] = np.nan
    with pytest.warns(UserWarning):
        empty = load_limb_sounders("2022-10-15", "ROMEX", file_path=tmp_path / "romex",
                                   original_z=True)
    assert "romex_testrmx_100m_v1a.nc" in empty.attrs["files_without_data"]
    # an Aeolus month file whose winds that day are all fill values
    folder = tmp_path / "aeolus"
    write_aeolus(folder)
    with netCDF4.Dataset(folder / "Aeolus_L2B_monthly_2019_01.nc", "a") as f:
        f["ray_hlos"][:] = 32767
    with pytest.warns(UserWarning):
        empty = load_limb_sounders("2019-01-10", "Aeolus", file_path=folder,
                                   original_z=True)
    assert "Aeolus_L2B_monthly_2019_01.nc" in empty.attrs["files_without_data"]
