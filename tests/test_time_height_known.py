"""Known-answer tests of load_time_height (ALOMAR reader), on small synthetic
HDF5 files laid out as the real ones (yearly files, one group per session,
named by its start date; [times x heights] arrays; POSIX times).

New code, no MATLAB equivalent. Every check is exact: values are only read,
placed and masked, never computed (altitudes are the lattice values rounded
to 0.01 km, compared exactly).

Sessions written (2007 file unless stated):

- A "2007-01-23": heights 30-105 km every 15 km, 4 rows at 13:00:30,
  13:05:30, 13:05:30 (repeated stamp), 13:10:30; a NaN inside row 0
  (missing), a NaN at the top of row 1 (no_observation), 50 K in row 2 at
  30 km and 450 K in row 3 at 90 km (removed), 450 K in row 3 at 105 km
  (kept: above 100 km);
- B "2007-01-24": heights 15-45 km, 4 rows every 2 h from 22:00:30, over
  midnight;
- C "2007-01-25": A's heights, 2 rows at 01:00:30 and 03:00:30, inside B
  (overlapping sessions);
- D "2007-12-31": A's heights, rows from 22:00:30 to 02:00:30 on 1 January;
- 2008 file: E "2008-01-05"; no 2009 file (nor 2006: the year before is
  looked in too, for sessions running into January).
"""

import warnings

import h5py
import numpy as np
import pytest

from atmosgw.curate.data import states, validate
from atmosgw.curate.obs import load_time_height
from atmosgw.curate.obs._limb_alomar import DATA_VERSION

Z_A = np.array([200, 300, 400, 500, 600, 700]) * 0.15  # 30 ... 105 km
Z_B = np.array([100, 200, 300]) * 0.15                 # 15, 30, 45 km


def posix(iso):
    return (np.datetime64(iso, "us") - np.datetime64("1970-01-01", "us")) / \
        np.timedelta64(1, "s")


def session(times, z, offset):
    """A session's arrays: distinct values T = offset + 10 * row + level."""
    t = np.array([posix(x) for x in times])
    n, m = t.size, z.size
    temp = offset + 10.0 * np.arange(n)[:, None] + np.arange(m)[None, :]
    unc = 1.0 + 0.01 * temp
    return {"time_lidar": t, "z_lidar": z, "T1_lidar": temp, "T1_lidar_unc": unc,
            "T2_lidar": temp + 0.5, "T2_lidar_unc": unc + 0.5,
            "T_ecmwf": temp - 1.0, "laser_power": 10.0 + np.arange(n),
            "odd": np.zeros((2, 2))}


def sessions_2007():
    a = session(["2007-01-23T13:00:30", "2007-01-23T13:05:30",
                 "2007-01-23T13:05:30", "2007-01-23T13:10:30"], Z_A, 200)
    a["T1_lidar"][0, 2] = np.nan   # inside row 0: missing
    a["T1_lidar"][1, 5] = np.nan   # top of row 1: no_observation
    a["T1_lidar"][2, 0] = 50.0     # removed
    a["T1_lidar"][3, 4] = 450.0    # 90 km: removed
    a["T1_lidar"][3, 5] = 450.0    # 105 km: kept (only checked below 100 km)
    for r, c in ((0, 2), (1, 5)):
        a["T1_lidar_unc"][r, c] = np.nan
    b = session(["2007-01-24T22:00:30", "2007-01-25T00:00:30",
                 "2007-01-25T02:00:30", "2007-01-25T04:00:30"], Z_B, 230)
    c = session(["2007-01-25T01:00:30", "2007-01-25T03:00:30"], Z_A, 260)
    d = session(["2007-12-31T22:00:30", "2008-01-01T00:00:30",
                 "2008-01-01T02:00:30"], Z_A, 280)
    return {"2007-01-23": a, "2007-01-24": b, "2007-01-25": c, "2007-12-31": d}


def write_year(folder, year, groups):
    with h5py.File(folder / f"{year}_{DATA_VERSION}_file4juli.h5", "w") as f:
        for name, arrays in groups.items():
            g = f.create_group(name)
            for key, value in arrays.items():
                g[key] = value


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    folder = tmp_path_factory.mktemp("alomar")
    write_year(folder, 2007, sessions_2007())
    write_year(folder, 2008, {"2008-01-05": session(["2008-01-05T12:00:30"], Z_A,
                                                    290)})
    return folder


def load(folder, time_range, **options):
    return load_time_height(time_range, "ALOMAR", file_path=folder, **options)


def flag_of(ds, state):
    return states.has_state(ds.temperature_flag.values, state)


# ------------------------------------------------------------------ values


def test_values_on_their_own_heights(files):
    ds = load(files, "2007-01-23", keep_outliers=True)
    validate(ds, "time_height")
    a = sessions_2007()["2007-01-23"]
    np.testing.assert_array_equal(ds.altitude.values, [30, 45, 60, 75, 90, 105])
    np.testing.assert_array_equal(ds.temperature.values, a["T1_lidar"])
    np.testing.assert_array_equal(ds.temperature_uncertainty.values,
                                  a["T1_lidar_unc"])
    assert ds.temperature.dtype == np.float64
    assert float(ds.latitude) == 69.38 and float(ds.longitude) == 16.08
    expected = np.array(["2007-01-23T13:00:30", "2007-01-23T13:05:30",
                         "2007-01-23T13:05:30", "2007-01-23T13:10:30"],
                        dtype="datetime64[us]")
    np.testing.assert_array_equal(ds.time.values, expected)


def test_union_grid_sorted_times_and_sessions(files):
    # B and C (overlapping sessions); A (2007-01-23) is read but has no rows
    ds = load(files, ["2007-01-24", "2007-01-25"], file_source=True)
    validate(ds, "time_height")
    np.testing.assert_array_equal(ds.altitude.values, [15, 30, 45, 60, 75, 90, 105])
    hours = (ds.time.values - np.datetime64("2007-01-24T00:00:30", "us")) \
        / np.timedelta64(1, "h")
    np.testing.assert_array_equal(hours, [22, 24, 25, 26, 27, 28])
    np.testing.assert_array_equal(ds.session.values, [1, 1, 2, 1, 2, 1])
    # rows counted through the file: A 1-4, B 5-8, C 9-10
    np.testing.assert_array_equal(ds.source_profile.values, [5, 6, 9, 7, 10, 8])
    np.testing.assert_array_equal(ds.source_file.values, [1] * 6)
    s = sessions_2007()
    b, c = s["2007-01-24"]["T1_lidar"], s["2007-01-25"]["T1_lidar"]
    np.testing.assert_array_equal(ds.temperature.values[[0, 1, 3, 5], :3], b)
    np.testing.assert_array_equal(ds.temperature.values[[2, 4], 1:], c)
    # heights a session does not have: no_observation, NaN
    lacks = np.zeros((6, 7), dtype=bool)
    lacks[[0, 1, 3, 5], 3:] = True
    lacks[[2, 4], 0] = True
    np.testing.assert_array_equal(flag_of(ds, "no_observation"), lacks)
    assert np.all(np.isnan(ds.temperature.values[lacks]))
    assert ds.attrs["state_counts"] == "temperature: no_observation 18"


# ------------------------------------------------------------------ states


def test_states(files):
    ds = load(files, "2007-01-23")
    t = ds.temperature.values
    assert flag_of(ds, "missing").sum() == 1 and flag_of(ds, "missing")[0, 2]
    assert flag_of(ds, "no_observation").sum() == 1
    assert flag_of(ds, "no_observation")[1, 5]
    bad = flag_of(ds, "outside_physical_range")
    assert bad.sum() == 2 and bad[2, 0] and bad[3, 4]
    assert np.all(np.isnan(t[bad]))
    assert np.all(np.isnan(ds.temperature_uncertainty.values[bad]))
    assert t[3, 5] == 450.0  # above 100 km: not checked against 400 K
    np.testing.assert_array_equal(np.isnan(ds.temperature_uncertainty.values),
                                  np.isnan(t))
    assert ds.attrs["state_counts"] == ("temperature: no_observation 1, missing 1, "
                                        "outside_physical_range 2")
    kept = load(files, "2007-01-23", keep_outliers=True)
    assert np.all(flag_of(kept, "outside_physical_range") == bad)
    assert kept.temperature.values[2, 0] == 50.0
    assert kept.temperature.values[3, 4] == 450.0
    assert np.isfinite(kept.temperature_uncertainty.values[2, 0])


# ------------------------------------------------------------------ options


def test_time_handling(files):
    strict = load(files, ["2007-01-23T13:00:00", "2007-01-23T13:06:00"],
                  time_handling="strict")
    assert strict.sizes["time"] == 3  # the repeated stamp included
    with pytest.warns(UserWarning, match="exactly at that instant"):
        none = load(files, "2007-01-23T13:01:00", time_handling="strict")
    assert none.sizes["time"] == 0
    generous = load(files, "2007-01-25T02:00:00", time_handling="generous")
    assert generous.sizes["time"] == 5  # all of 25 January: B 3 rows, C 2


def test_height_range(files):
    ds = load(files, "2007-01-23", height_range=(40, 95))
    np.testing.assert_array_equal(ds.altitude.values, [45, 60, 75, 90])
    assert ds.temperature.shape == (4, 4)


def test_additional_vars(files):
    ds = load(files, "2007-01-23", additional_vars=["T2_lidar", "laser_power"])
    a = sessions_2007()["2007-01-23"]
    np.testing.assert_array_equal(ds.T2_lidar.values, a["T2_lidar"])
    assert ds.T2_lidar.dims == ("time", "level")
    np.testing.assert_array_equal(ds.laser_power.values, a["laser_power"])
    assert ds.laser_power.dims == ("time",)
    with pytest.raises(KeyError, match="Nope"):
        load(files, "2007-01-23", additional_vars=["Nope"])
    with pytest.raises(ValueError, match="only"):
        load(files, "2007-01-23", additional_vars=["odd"])
    with pytest.raises(ValueError, match="not an additional variable"):
        load(files, "2007-01-23", additional_vars=["z_lidar"])


def test_variable_per_time_in_one_session_per_point_in_another(tmp_path):
    # (independent review, 2026-11-26: was broadcast over levels, no error)
    a = session(["2007-02-01T12:00:30", "2007-02-01T12:05:30"], Z_A, 200)
    b = session(["2007-02-01T13:00:30"], Z_A, 200)
    a["laser_power"] = np.ones((2, Z_A.size))
    write_year(tmp_path, 2007, {"2007-02-01a": a, "2007-02-01b": b})
    with pytest.raises(ValueError, match="per time in some sessions"):
        load(tmp_path, "2007-02-01", additional_vars=["laser_power"])


def test_look_back_into_previous_year(files):
    # D started on 2007-12-31, in the 2007 file; its 1 January rows belong to
    # a request for 2008-01-01
    ds = load(files, "2008-01-01", file_source=True)
    assert ds.sizes["time"] == 2
    assert list(ds.original_files.values) == [
        f"2007_{DATA_VERSION}_file4juli.h5", f"2008_{DATA_VERSION}_file4juli.h5"]
    np.testing.assert_array_equal(ds.source_file.values, [1, 1])
    np.testing.assert_array_equal(ds.source_profile.values, [12, 13])
    assert ds.attrs["files_not_found"] == ""


def test_year_without_file_and_empty_output(files):
    ds = load(files, "2009-03-01")
    validate(ds, "time_height")
    assert ds.sizes["time"] == 0 and ds.sizes["level"] == 0
    assert ds.attrs["files_not_found"] == "2009"
    assert ds.attrs["files_read"] == ""


def test_options_refused_and_warnings(files):
    with pytest.raises(TypeError, match="load_limb_sounders"):
        load(files, "2007-01-23", original_z=True)
    with pytest.raises(TypeError, match="unexpected"):
        load(files, "2007-01-23", nonsense=1)
    with pytest.raises(ValueError, match="Unknown instrument"):
        load_time_height("2007-01-23", "SABER", file_path=files)
    with pytest.raises(TypeError, match="True or False"):
        load(files, "2007-01-23", keep_outliers=1)
    with pytest.warns(UserWarning, match="only available"):
        load(files, "2001-01-01")
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # nothing to warn about in normal use
        load(files, "2007-01-23")


def test_off_lattice_heights_refused(tmp_path):
    bad = session(["2007-02-01T12:00:30"], np.array([30.0, 30.07]), 200)
    write_year(tmp_path, 2007, {"2007-02-01": bad})
    with pytest.raises(ValueError, match="0.15 km grid"):
        load(tmp_path, "2007-02-01")


def test_verbose_and_attributes(files, capsys):
    ds = load(files, "2007-01-23", verbose=True)
    assert "Loading ALOMAR data: 2007-01-23" in capsys.readouterr().out
    assert ds.attrs["format"] == "time_height"
    assert ds.attrs["format_version"] == 1
    assert ds.attrs["site_name"] == "ALOMAR, Andoya, Norway"
    assert ds.attrs["heights"] == "original"
    assert ds.attrs["data"] == "raw" and ds.attrs["preprocessing"] == ""
    assert ds.attrs["source"] == "atmosgw.curate.obs.load_time_height"
    assert ds.attrs["availability"] == (
        "eastward_wind: not_in_files; northward_wind: not_in_files; "
        "upward_wind: not_in_files; pressure: not_in_files; "
        "source_profile: not_requested; source_file: not_requested; "
        "original_files: not_requested")
    assert ds.temperature.attrs["source_variable"] == "T1_lidar"
    assert "UNVERIFIED" in ds.temperature_uncertainty.attrs["comment"]


def test_long_session_found_by_its_times(tmp_path):
    # sessions can last days (8.05 in the real files): rows on 8 February come
    # from the session named 2007-02-01, and one named in December of the year
    # before can run into January
    long = session(["2007-02-01T12:00:30", "2007-02-05T12:00:30",
                    "2007-02-08T12:00:30"], Z_A, 200)
    write_year(tmp_path, 2007, {"2007-02-01": long})
    write_year(tmp_path, 2006, {"2006-12-30": session(
        ["2006-12-30T20:00:30", "2007-01-06T01:00:30"], Z_B, 210)})
    ds = load(tmp_path, "2007-02-08", file_source=True)
    np.testing.assert_array_equal(ds.source_profile.values, [3])
    ds = load(tmp_path, "2007-01-06", file_source=True)
    np.testing.assert_array_equal(ds.source_profile.values, [2])
    assert list(ds.original_files.values) == [f"2006_{DATA_VERSION}_file4juli.h5",
                                              f"2007_{DATA_VERSION}_file4juli.h5"]
