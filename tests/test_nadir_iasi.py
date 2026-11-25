"""Known-answer tests for the IASI reader and screen (synthetic files).

Synthetic orbit files are written to a temporary folder in the layout the
reader expects; no real data are needed. Tolerances: exact (states, row
counts, values that are passed through unchanged).
"""

import netCDF4
import numpy as np

from atmosgw.curate.data import states, validate
from atmosgw.curate.obs import load_nadir_sounders
from atmosgw.curate.obs._nadir_iasi import granule_of, screen

DAY_SECONDS = ((np.datetime64("2007-11-07") - np.datetime64("2000-01-01"))
               / np.timedelta64(1, "s"))
Z = np.r_[0:61:3, 65:91:5].astype(np.float32)


def write_orbit(root, start_seconds, n_rows, temp=None, with_temp=True):
    """One IASI orbit file starting at start_seconds (rows every 4 s, 60 across)."""
    t = np.datetime64("2000-01-01T00:00:00") + np.timedelta64(int(start_seconds), "s")
    stamp = str(t).replace("-", "").replace(":", "").replace("T", "")
    day = str(t)[:10]
    doy = int((np.datetime64(day) - np.datetime64(day[:4] + "-01-01")) /
              np.timedelta64(1, "D")) + 1
    folder = root / day[:4] / f"{doy:03d}"
    folder.mkdir(parents=True, exist_ok=True)
    rows = np.arange(n_rows)
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


def test_screen_function_matches_definition():
    x = np.arange(1.0, 101.0)[:, None, None] * np.ones((1, 1, 2))
    x[0, 0, 1] = 1000.0
    out = screen(x)
    assert out[0, 0, 1] and out.sum() == 1
