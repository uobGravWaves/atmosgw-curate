"""Known-answer tests for the IASI reader and screen (synthetic files).

Synthetic orbit files are written to a temporary folder in the layout the
reader expects; no real data are needed. Tolerances: exact (states, row
counts, values that are passed through unchanged).
"""

import netCDF4
import numpy as np

from atmosgw.curate.data import states, validate
from atmosgw.curate.obs import load_nadir_sounders
from atmosgw.curate.obs._nadir_iasi import granule_of
from atmosgw.curate.obs._screen import pooled_screen, resolve, rolling_screen

DAY_SECONDS = ((np.datetime64("2007-11-07") - np.datetime64("2000-01-01"))
               / np.timedelta64(1, "s"))
Z = np.r_[0:61:3, 65:91:5].astype(np.float32)


def write_orbit(root, start_seconds, n_rows, temp=None, with_temp=True, skip=None):
    """One IASI orbit file starting at start_seconds (rows every 4 s, 60
    across); skip (row, n): n rows absent before that row (times and
    positions advanced past them)."""
    t = np.datetime64("2000-01-01T00:00:00") + np.timedelta64(int(start_seconds), "s")
    stamp = str(t).replace("-", "").replace(":", "").replace("T", "")
    day = str(t)[:10]
    doy = int((np.datetime64(day) - np.datetime64(day[:4] + "-01-01")) /
              np.timedelta64(1, "D")) + 1
    folder = root / day[:4] / f"{doy:03d}"
    folder.mkdir(parents=True, exist_ok=True)
    rows = np.arange(n_rows)
    if skip is not None:
        rows = np.where(rows >= skip[0], rows + skip[1], rows)
    lat = np.broadcast_to(-40 + 0.05 * rows[:, None], (n_rows, 60)) + 0 * np.arange(60)
    lon = np.broadcast_to(20 + 0.3 * np.arange(60), (n_rows, 60)) + 0.01 * rows[:, None]
    time = start_seconds + np.broadcast_to(4.0 * rows[:, None], (n_rows, 60))
    if temp is None:
        temp = np.broadcast_to(220 + Z, (n_rows, 60, Z.size)).astype(np.float32)
    path = folder / f"iasi_{stamp}.nc"
    with netCDF4.Dataset(path, "w") as f:
        f.createDimension("L1_NTRACK", n_rows)
        f.createDimension("L1_NXTRACK", 60)
        f.createDimension("RET_NP", Z.size)
        for name, value in (("l1_time", time), ("l1_lat", lat), ("l1_lon", lon)):
            f.createVariable(name, "f8", ("L1_NTRACK", "L1_NXTRACK"))[:] = value
        f.createVariable("ret_z", "f4", ("RET_NP",))[:] = Z
        if with_temp:
            f.createVariable("ret_temp", "f4",
                             ("L1_NTRACK", "L1_NXTRACK", "RET_NP"))[:] = temp
    return path


def test_granule_bins_exact_at_edges():
    t = DAY_SECONDS + np.array([0, 359.999, 360, 86399.999])
    day, g = granule_of(t)
    assert g.tolist() == [1, 1, 2, 240]
    assert len(set(day.tolist())) == 1


def test_granule_spanning_a_file_start_gets_both_files(tmp_path):
    # file 1: 00:00 - 00:50:56; file 2 from 00:50:56; granule 9 is 00:48 - 00:54
    write_orbit(tmp_path, DAY_SECONDS, int(3056 / 4))
    write_orbit(tmp_path, DAY_SECONDS + 3056, 900)
    ds = load_nadir_sounders(None, "IASI", granules=[("2007-11-07", 9)],
                             file_path=tmp_path, original_grid=True)
    assert ds.sizes["along_track"] == 90  # 6 min of rows every 4 s
    assert ds.attrs["files_read"].count("iasi_") == 2
    assert set(ds.granule.values) == {9}


def test_short_dropout_splits_the_swath(tmp_path):
    # 5 rows (20 s) missing after the first 40: with the 12 s segment rule the
    # regular grid has two segments and no value dated inside the gap (the
    # earlier 180 s rule bridged it with good-looking values; found in review)
    write_orbit(tmp_path, DAY_SECONDS, 180, skip=(40, 5))
    ds = load_nadir_sounders(None, "IASI", granules=[("2007-11-07", 1),
                                                     ("2007-11-07", 2)],
                             file_path=tmp_path, screen="none")
    validate(ds, "swath")
    assert len(np.unique(ds.segment.values)) == 2
    sec = (ds.time.values - np.datetime64("2007-11-07")) / np.timedelta64(1, "s")
    finite = np.isfinite(ds.temperature.values[..., 0])
    assert not np.any(finite & (sec > 4 * 39 + 2) & (sec < 4 * 45 - 2))


def test_stray_pixel_time_or_short_step_does_not_split(tmp_path):
    # one pixel 590 s off its scan (kept: within the 600 s stray rule) moved
    # a row's mean time by 9.8 s, enough to fake a 12 s gap; the row's
    # median ignores it. A step of 11 s (under 1.5 scans) is not a gap.
    path = write_orbit(tmp_path, DAY_SECONDS, 90)
    with netCDF4.Dataset(path, "a") as f:
        t = f["l1_time"][:]
        t[50, 7] += 590.0
        t[60:] += 7.0  # (rows 59 -> 60: 11 s)
        f["l1_time"][:] = t
    ds = load_nadir_sounders(None, "IASI", granules=[("2007-11-07", 1)],
                             file_path=tmp_path, screen="none")
    assert len(np.unique(ds.segment.values)) == 1


def test_file_without_temperature_is_skipped_and_recorded(tmp_path):
    write_orbit(tmp_path, DAY_SECONDS, 900, with_temp=False)
    ds = load_nadir_sounders(None, "IASI", granules=[("2007-11-07", 3)],
                             file_path=tmp_path, original_grid=True)
    assert ds.sizes["along_track"] == 0
    assert "no 'ret_temp'" in ds.attrs["files_not_found"]
    validate(ds, "swath")


def test_screen_removes_spike_keeps_wave(tmp_path):
    rng = np.random.default_rng(1)
    rows = np.arange(900)
    wave = 6 * np.cos(2 * np.pi * (rows - 135) / 12) * np.exp(-((rows - 135) / 6) ** 2)
    temp = (220 + Z[None, None, :] + 0.25 * np.arange(60)[None, :, None]
            + 1.5 * rng.standard_normal((900, 60, Z.size))
            + wave[:, None, None] * (np.arange(60) < 30)[None, :, None])
    temp = temp.astype(np.float32)
    temp[200, 45, 13] += 100
    write_orbit(tmp_path, DAY_SECONDS, 900, temp)
    wanted = [("2007-11-07", k) for k in (1, 2, 3)]
    ds = load_nadir_sounders(None, "IASI", granules=wanted,
                             file_path=tmp_path, original_grid=True)
    flag = ds.temperature_flag.values
    outlier = states.has_state(flag, "statistical_outlier")
    assert outlier.sum() == 1 and outlier[200, 45, 13]
    assert np.isnan(ds.temperature.values[200, 45, 13])
    others = ~outlier
    np.testing.assert_array_equal(ds.temperature.values[others], temp[:270][others])


def test_pooled_screen_matches_definition():
    x = np.arange(1.0, 101.0)[:, None, None] * np.ones((1, 1, 2))
    x[0, 0, 1] = 1000.0
    out = pooled_screen(x, None, None, None)
    assert out[0, 0, 1] and out.sum() == 1


# ------------------------------------------------- rolling screen, swappable


def row_times(n_rows, start=733353.0, seconds=4.0, n_cross=60):
    """[cross x along] datenums, one row every `seconds`."""
    rows = np.arange(n_rows) * seconds / 86400
    return start + np.broadcast_to(rows, (n_cross, n_rows))


def test_rolling_screen_follows_slow_changes():
    """A 30 K drift over an hour is not flagged; a 100 K spike is."""
    n = 900
    rng = np.random.default_rng(2)
    drift = np.linspace(0, 30, n)
    temp = (220 + drift[None, :, None] + rng.standard_normal((60, n, 3))).copy()
    temp[10, 450, 1] += 100
    lat = np.zeros((60, n))
    out = rolling_screen(temp, row_times(n), lat, lat, split_day_night=False)
    assert out.sum() == 1 and out[10, 450, 1]
    # the pooled (MATLAB) screen over the hour would see the drift as spread
    assert pooled_screen(temp, None, None, None)[10, 450, 1]


def test_rolling_screen_day_night_groups():
    """A 15 K day/night offset with a thin night sliver: split keeps the sliver."""
    n = 300
    rng = np.random.default_rng(3)
    temp = 220 + rng.standard_normal((60, n, 1))
    time = row_times(n, start=733353.5)  # noon UTC
    lon = np.broadcast_to(np.linspace(-30, 30, 60)[:, None], (60, n)).copy()
    lat = np.zeros((60, n))
    lon[:6] = 180.0  # midnight longitude: night points
    temp[:6] += 15
    split = rolling_screen(temp, time, lat, lon, min_points=50)
    pooled = rolling_screen(temp, time, lat, lon, split_day_night=False)
    assert split.sum() == 0
    assert pooled[:6].mean() > 0.9  # the whole sliver flagged without the split


def test_rolling_screen_too_few_points_not_screened():
    temp = 220 + np.random.default_rng(5).standard_normal((60, 10, 1))
    temp[0, 5, 0] = 500
    out = rolling_screen(temp, row_times(10), np.zeros((60, 10)), np.zeros((60, 10)),
                         split_day_night=False)  # 600 points >= 500: screened
    assert out[0, 5, 0]
    out = rolling_screen(temp, row_times(10), np.zeros((60, 10)), np.zeros((60, 10)),
                         split_day_night=False, min_points=1000)
    assert not out.any()


def test_flags_do_not_depend_on_the_request(tmp_path):
    rng = np.random.default_rng(4)
    temp = (220 + Z[None, None, :] + 1.5 * rng.standard_normal((900, 60, Z.size)))
    temp[100:110, 5, 10] += 40  # moderate outliers near a granule edge
    write_orbit(tmp_path, DAY_SECONDS, 900, temp.astype(np.float32))
    one = load_nadir_sounders(None, "IASI", granules=[("2007-11-07", 2)],
                              file_path=tmp_path, original_grid=True)
    three = load_nadir_sounders(None, "IASI",
                                granules=[("2007-11-07", k) for k in (1, 2, 3)],
                                file_path=tmp_path, original_grid=True)
    middle = three.granule.values == 2
    np.testing.assert_array_equal(one.temperature_flag.values,
                                  three.temperature_flag.values[middle])
    assert one.attrs["files_not_found"] == ""  # padding granules are not requests


def test_exactly_100k_is_outside_physical_range(tmp_path):
    # IASI's retrieval floor: removed whatever the screen; 100.001 K is kept
    temp = np.broadcast_to(220 + Z, (270, 60, Z.size)).astype(np.float32).copy()
    temp[50, 10, :] = 100.0           # a whole profile at the floor
    temp[70, 20, 12] = np.float32(100.001)
    write_orbit(tmp_path, DAY_SECONDS, 270, temp)
    ds = load_nadir_sounders(None, "IASI", granules=[("2007-11-07", 1)],
                             file_path=tmp_path, original_grid=True, screen="none")
    bad = states.has_state(ds.temperature_flag.values, "outside_physical_range")
    assert bad[50, 10, :].all() and bad.sum() == Z.size
    assert np.all(np.isnan(ds.temperature.values[50, 10, :]))
    assert ds.temperature.values[70, 20, 12] == np.float32(100.001)


def test_screen_is_swappable(tmp_path):
    write_orbit(tmp_path, DAY_SECONDS, 900)
    args = dict(granules=[("2007-11-07", 2)], file_path=tmp_path, original_grid=True)

    def bottom_level(temp, time, lat, lon):
        out = np.zeros(np.shape(temp), dtype=bool)
        out[..., 0] = True
        return out

    custom = load_nadir_sounders(None, "IASI", screen=bottom_level, **args)
    outlier = states.has_state(custom.temperature_flag.values, "statistical_outlier")
    assert outlier[..., 0].all() and not outlier[..., 1:].any()
    assert custom.attrs["screen"] == "custom: bottom_level"
    none = load_nadir_sounders(None, "IASI", screen="none", **args)
    flags = none.temperature_flag.values
    assert not states.has_state(flags, "statistical_outlier").any()
    assert none.attrs["screen"] == "none"
    tuned = load_nadir_sounders(None, "IASI", screen=("rolling", {"factor": 5}), **args)
    assert tuned.attrs["screen"] == "rolling(factor=5)"
    default = load_nadir_sounders(None, "IASI", **args)
    assert default.attrs["screen"] == "rolling()"
    assert resolve(None, ("none", {}))[0] is None
    import pytest
    with pytest.raises(ValueError):
        load_nadir_sounders(None, "IASI", screen="median", **args)
