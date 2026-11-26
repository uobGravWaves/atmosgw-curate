"""Known-answer tests for fill values and bad times in nadir files (AIRS, IASI).

The files declare no fill values; obs._common.nadir_fill_to_nan sets them to
NaN (census of the eepc-0184 archive, 2026-11-26): -999 (padded AIRS daily-file
granules), -9999 AIRS positions, latitude = longitude = 0 and zero times
(IASI, one AIRS granule), stray pixel times and misdated scans. Synthetic files
only; states, counts and NaN patterns are exact.
"""

import netCDF4
import numpy as np
import pytest
from test_nadir_iasi import DAY_SECONDS, write_orbit
from test_nadir_known import D, seconds, write_granule

from atmosgw.curate.data import states, validate
from atmosgw.curate.obs import load_nadir_sounders
from atmosgw.curate.obs._common import nadir_fill_to_nan


def edit(path, **changes):
    """Overwrite parts of a synthetic file: name=(index, value)."""
    with netCDF4.Dataset(path, "a") as f:
        f.set_auto_maskandscale(False)
        for name, (index, value) in changes.items():
            v = f.variables[name][...]
            v[index] = value
            f.variables[name][...] = v


def airs_path(root, g):
    return root / "2007" / "001" / f"airs_2007_001_{g:03d}.nc"


def test_fill_rules():
    lat = np.full((4, 5), 10.0)
    lon = np.full((4, 5), 20.0)
    time = 1e8 + np.arange(4)[:, None] * 8 + np.zeros((1, 5))
    temp = np.full((4, 5, 2), 250.0)
    lat[0, 0], lon[0, 1] = -999, -9999   # fill positions (both lat and lon go)
    lat[1, 2], lon[1, 2] = 0, 0          # 0/0: no geolocation
    lat[1, 3] = 0                        # latitude 0 alone is real
    time[2, 4] = 0                       # zero time
    time[2, 0] = 1e8 + 20 * 86400        # stray pixel
    time[3, :] += 5 * 86400              # misdated scan
    temp[0, 0, 0], temp[0, 1, 1] = -999, -9999
    counts = nadir_fill_to_nan(lat, lon, time, temp, (1e8, 1e8 + 86400))
    assert counts == {"position": 3, "time": 1, "stray": 1, "misdated": 1, "temp": 2}
    assert np.isnan(lat[0, 0]) and np.isnan(lon[0, 0]) and np.isnan(lat[0, 1])
    assert np.isnan(lat[1, 2]) and np.isnan(lon[1, 2]) and lat[1, 3] == 0
    assert np.isnan(time[2, 4]) and np.isnan(time[2, 0])
    assert np.all(np.isnan(time[3])) and np.all(np.isfinite(time[:2]))
    assert np.isnan(temp).sum() == 2


def test_airs_fill_granule_recorded(tmp_path):
    for g in (53, 54, 55):
        write_granule(tmp_path, g, seconds(g))
    edit(airs_path(tmp_path, 54), l1_lat=(..., -999), l1_lon=(..., -999),
         l1_time=(..., -999), ret_temp=(..., -999))
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, g) for g in (53, 54, 55)],
                             file_path=tmp_path)
    validate(ds, "swath")
    assert sorted(set(ds.granule.values)) == [53, 55]
    assert "2007-01-01 granule 54 (no valid times)" in ds.attrs["files_not_found"]
    # a time range inside the fill granule: empty, but recorded (nominal time)
    t = np.datetime64(D) + np.timedelta64(int(seconds(54) - seconds(1)) + 328, "s")
    span = [t + np.timedelta64(60, "s"), t + np.timedelta64(300, "s")]
    ds = load_nadir_sounders(span, "AIRS", file_path=tmp_path)
    assert ds.sizes["along_track"] == 0
    assert "granule 54 (no valid times)" in ds.attrs["files_not_found"]


def test_airs_missing_scans_not_regridded(tmp_path):
    for g in (53, 54):
        write_granule(tmp_path, g, seconds(g))
    # granule 54: scans 1-60 without positions (-9999), times kept
    edit(airs_path(tmp_path, 54), l1_lat=(slice(0, 60), -9999),
         l1_lon=(slice(0, 60), -9999))
    native = load_nadir_sounders(None, "AIRS", granules=[(D, 53), (D, 54)],
                                 file_path=tmp_path, original_grid=True)
    validate(native, "swath")
    assert native.sizes["along_track"] == 270
    bad = states.has_state(native.temperature_flag.values, "bad_geolocation")
    rows = np.flatnonzero(bad.all(axis=(1, 2)))
    assert rows.tolist() == list(range(135, 195))  # rows 1-60 of granule 54
    assert not bad[:135].any() and not bad[195:].any()
    assert sorted(set(native.segment.values)) == [1, 2, 3]
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53), (D, 54)],
                             file_path=tmp_path)
    validate(ds, "swath")
    assert sorted(set(ds.segment.values)) == [1, 3]  # segment 2 not regridded
    assert np.all(np.isfinite(ds.latitude.values))


def test_airs_stray_time_is_bad_geolocation(tmp_path):
    write_granule(tmp_path, 53, seconds(53))
    edit(airs_path(tmp_path, 53), l1_time=((70, 10), seconds(53) + 20 * 86400))
    native = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                                 original_grid=True)
    bad = states.has_state(native.temperature_flag.values, "bad_geolocation")
    assert bad.any(axis=2).sum() == 1 and bad[70, 10].all()
    assert np.isnat(native.time.values[70, 10])
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path)
    validate(ds, "swath")
    assert sorted(set(ds.segment.values)) == [1]  # a stray time alone: regridded


def test_padding_neighbour_fill_not_reported(tmp_path):
    for g in (52, 53, 54):
        write_granule(tmp_path, g, seconds(g))
    edit(airs_path(tmp_path, 54), l1_lat=(..., -999), l1_lon=(..., -999),
         l1_time=(..., -999), ret_temp=(..., -999))
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                             screen="rolling")
    assert "granule 54" not in ds.attrs["files_not_found"]
    assert sorted(set(ds.granule.values)) == [53]


def test_no_placeable_scans_warns_native(tmp_path):
    write_granule(tmp_path, 53, seconds(53))
    edit(airs_path(tmp_path, 53), l1_lat=(..., -9999), l1_lon=(..., -9999))
    with pytest.warns(UserWarning, match="native grid is returned"):
        ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path)
    assert ds.attrs["grid"] == "native" and ds.sizes["along_track"] == 135


def test_iasi_zero_scans_kept_in_their_granule(tmp_path):
    path = write_orbit(tmp_path, DAY_SECONDS, 900)  # 00:00-01:00, rows every 4 s
    # rows 100-149 (granule 2: 00:06:40-00:09:56): no geolocation, zero times at
    # the swath edges (as in the files); rows 300-309: every time zero
    edit(path, l1_lat=(slice(100, 150), 0), l1_lon=(slice(100, 150), 0))
    edit(path, l1_time=((slice(100, 150), [0, 1, 2, 3, 58, 59]), 0))
    edit(path, l1_time=(slice(300, 310), 0))
    native = load_nadir_sounders(None, "IASI", granules=[("2007-11-07", 2)],
                                 file_path=tmp_path, original_grid=True)
    assert native.sizes["along_track"] == 90  # none lost to a 2000-01-01 granule
    bad = states.has_state(native.temperature_flag.values, "bad_geolocation")
    assert bad.all(axis=(1, 2)).sum() == 50
    assert native.attrs["files_not_found"] == ""  # granule 4's loss: not asked
    ds = load_nadir_sounders(None, "IASI", granules=[("2007-11-07", 2)],
                             file_path=tmp_path)
    validate(ds, "swath")
    assert np.all(np.isfinite(ds.latitude.values))
    assert len(set(ds.segment.values)) == 2  # rows 100-149 left out between them
    # scans with no valid time at all: dropped and recorded
    g4 = load_nadir_sounders(None, "IASI", granules=[("2007-11-07", 4)],
                             file_path=tmp_path, original_grid=True)
    assert g4.sizes["along_track"] == 80
    assert "2007-11-07 granule 4: 10 scans without valid times" in (
        g4.attrs["files_not_found"])
