"""Known-answer and format tests for load_nadir_sounders (AIRS).

Synthetic granule files (written to a temporary folder) need no real data:
data states, outlier handling, missing granules and the swath format.
Real-data checks (AIRS 2007-01-01/02) are skipped if those files are absent.

Tolerances: positions and distances 1e-6 km or degrees where stated; states
and counts exact.
"""

import inspect
import warnings

import netCDF4
import numpy as np
import pytest

from atmosgw.curate.data import problems, states, validate
from atmosgw.curate.geodesy import km2deg, reckon
from atmosgw.curate.obs import load_nadir_sounders
from atmosgw.curate.paths import DataDirNotSetError, dataset_dir

D = "2007-01-01"


# ------------------------------------------------------------ synthetic files


def write_granule(root, granule, start_seconds, temp=None, cross_scale=1.0):
    """An AIRS granule file (135 along x 90 cross x 27 levels) on a curved track
    (cross_scale: scale of the cross-track offsets)."""
    folder = root / "2007" / "001"
    folder.mkdir(parents=True, exist_ok=True)
    n_al, n_x, n_z = 135, 90, 27
    head = 350 + np.linspace(0, 4, n_al)
    mid_lat, mid_lon = np.empty(n_al), np.empty(n_al)
    mid_lat[0], mid_lon[0] = 10 + 18 * 135 * (granule - 1) / 111.2, 30.0
    for k in range(1, n_al):
        mid_lat[k], mid_lon[k] = reckon(mid_lat[k - 1], mid_lon[k - 1], km2deg(18),
                                        head[k - 1])
    off = (cross_scale * np.sign(np.arange(n_x) - 44.5) * 13.5
           * np.abs(np.arange(n_x) - 44.5) ** 1.1)
    lat, lon = reckon(np.broadcast_to(mid_lat[:, None], (n_al, n_x)),
                      np.broadcast_to(mid_lon[:, None], (n_al, n_x)),
                      km2deg(np.broadcast_to(off, (n_al, n_x))),
                      np.broadcast_to(head[:, None] + 90, (n_al, n_x)))
    z = np.r_[0:61:3, 65:91:5].astype(np.float32)
    if temp is None:
        temp = np.broadcast_to(220 + 0.5 * z, (n_al, n_x, n_z)).astype(np.float32)
    time = start_seconds + np.broadcast_to(np.arange(n_al)[:, None] * 8 / 3,
                                           (n_al, n_x))
    path = folder / f"airs_2007_001_{granule:03d}.nc"
    with netCDF4.Dataset(path, "w") as f:
        f.createDimension("L1_NTRACK", n_al)
        f.createDimension("L1_NXTRACK", n_x)
        f.createDimension("RET_NP", n_z)
        for name, value in (("l1_lat", lat), ("l1_lon", lon), ("l1_time", time)):
            f.createVariable(name, "f8", ("L1_NTRACK", "L1_NXTRACK"))[:] = value
        v = f.createVariable("ret_temp", "f4", ("L1_NTRACK", "L1_NXTRACK", "RET_NP"),
                             fill_value=np.float32(-999))
        v[:] = temp
        f.createVariable("ret_z", "f4", ("RET_NP",))[:] = z
    return temp


def seconds(granule):
    """Nominal start of a 2007-01-01 granule, seconds since 2000-01-01."""
    day = (np.datetime64(D) - np.datetime64("2000-01-01")) / np.timedelta64(1, "s")
    return day + 5 * 60 + 28 + (granule - 1) * 360


def test_states_outliers_and_format(tmp_path):
    temp = np.broadcast_to(220 + 0.5 * np.r_[0:61:3, 65:91:5],
                           (135, 90, 27)).astype(np.float32).copy()
    temp[60, 45, 8] = 50.0     # outside the physical range
    temp[100, 20, 4] = -999.0  # fill value: missing
    write_granule(tmp_path, 53, seconds(53), temp)
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                             original_z=True)
    validate(ds, "swath")
    flag = ds.temperature_flag.values
    assert states.has_state(flag, "outside_physical_range")[..., 8].sum() > 0
    assert states.has_state(flag, "missing")[..., 4].sum() > 0
    assert np.all(np.isnan(ds.temperature.values[states.no_value_mask(flag)]))
    # other levels untouched, and exactly the file values (a linear field)
    level0 = ds.temperature.values[..., 0]
    assert np.all(level0 == np.float32(220))
    # keep_outliers: the 50 K point is used (values near it are below 220 + 4)
    kept = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                               original_z=True, keep_outliers=True)
    near = states.has_state(kept.temperature_flag.values, "outside_physical_range")
    assert np.all(np.isfinite(kept.temperature.values[near]))
    assert np.nanmin(kept.temperature.values[..., 8]) < 220 + 0.5 * 24


def test_airs_keeps_exactly_100k(tmp_path):
    # only IASI treats 100 K as a retrieval floor; AIRS keeps it (as MATLAB)
    temp = np.broadcast_to(220 + 0.5 * np.r_[0:61:3, 65:91:5],
                           (135, 90, 27)).astype(np.float32).copy()
    temp[60, 45, 8] = 100.0
    write_granule(tmp_path, 53, seconds(53), temp)
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                             original_grid=True)
    assert ds.temperature.values[60, 45, 8] == 100
    flag = ds.temperature_flag.values
    assert not states.has_state(flag, "outside_physical_range").any()


def test_missing_granule_gives_two_segments(tmp_path):
    for g in (53, 55):
        write_granule(tmp_path, g, seconds(g))
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53), (D, 54), (D, 55)],
                             file_path=tmp_path)
    validate(ds, "swath")
    assert sorted(set(ds.segment.values)) == [1, 2]
    assert sorted(set(ds.granule.values)) == [53, 55]
    assert "2007-01-01 granule 54" in ds.attrs["files_not_found"]
    assert ds.sizes["along_track"] == 2 * 135
    # each segment's along-track distance starts at 0
    for s in (1, 2):
        assert ds.distance_along.values[ds.segment.values == s][0] == 0


def test_time_selection_by_mean_time(tmp_path):
    for g in (50, 51, 52, 53):
        write_granule(tmp_path, g, seconds(g))
    ds = load_nadir_sounders([f"{D}T05:00", f"{D}T05:20"], "AIRS", file_path=tmp_path)
    # granule g's mean time is 00:05:28 + 6(g - 1) min + 3 min: 50-52 lie in range
    assert sorted(set(ds.granule.values)) == [50, 51, 52]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(UserWarning, match="strict"):
            load_nadir_sounders(f"{D}T05:00", "AIRS", file_path=tmp_path,
                                time_handling="strict")


def test_optional_outputs_and_availability(tmp_path):
    write_granule(tmp_path, 53, seconds(53))
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                             file_source=True, day_night=True,
                             keep_native_geolocation=True)
    assert problems(ds, "swath") == []
    assert list(ds.original_files.values) == ["airs_2007_001_053.nc"]
    assert ds.native_latitude.shape == (135, 90)
    plain = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path)
    avail = states.parse_availability(plain.attrs["availability"])
    assert avail["satellite_latitude"] == "not_in_files"
    assert avail["source_file"] == "not_requested"
    assert avail["temperature_uncertainty"] == "not_provided"


def test_lat_range_keeps_whole_rows(tmp_path):
    write_granule(tmp_path, 53, seconds(53))
    full = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path)
    mid = float(np.median(full.latitude))
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                             lat_range=(mid, 90))
    assert 0 < ds.sizes["along_track"] < full.sizes["along_track"]
    assert ds.sizes["cross_track"] == 90
    assert np.all(np.max(ds.latitude.values, axis=1) >= mid)


def test_no_files(tmp_path):
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path)
    assert ds.sizes["along_track"] == 0
    assert "granule 53" in ds.attrs["files_not_found"]
    validate(ds, "swath")


def test_rejects():
    with pytest.raises(NotImplementedError):
        load_nadir_sounders(D, "IASI", product="2d", file_path="/nonexistent")
    with pytest.raises(ValueError):
        load_nadir_sounders(D, "AIRS", interp_type="natural", file_path="/nonexistent")
    with pytest.raises(ValueError):
        load_nadir_sounders(D, "IASI-NG", file_path="/nonexistent")
    with pytest.raises(ValueError):
        load_nadir_sounders(None, "AIRS", file_path="/nonexistent")


# ------------------------------------------------------------ real data


@pytest.fixture(scope="module")
def airs():
    try:
        folder = dataset_dir("AIRS")
    except (DataDirNotSetError, FileNotFoundError) as error:
        pytest.skip(f"no AIRS data: {error}")
    if not (folder / "2007" / "airs_2007d001.nc").is_file():
        pytest.skip("AIRS 2007-01-01 not found")
    return folder


def test_real_daily_files_are_preprocessed(airs):
    # the daily files are rounded to whole kelvin; granule files are not
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], original_grid=True)
    assert ds.attrs["data"] == "preprocessed"
    assert "rounded to whole kelvin" in ds.attrs["preprocessing"]


def test_real_midnight_is_one_segment(airs):
    ds = load_nadir_sounders(None, "AIRS", granules=[(D, 240), ("2007-01-02", 1)])
    validate(ds, "swath")
    assert set(ds.segment.values) == {1}
    assert list(dict.fromkeys(ds.granule.values)) == [240, 1]
    days = np.unique(ds.granule_day.values.astype("datetime64[D]")).astype(str)
    assert list(days) == ["2007-01-01", "2007-01-02"]


def test_real_signature_matches_documented_defaults():
    defaults = {k: p.default for k, p in
                inspect.signature(load_nadir_sounders).parameters.items()
                if p.default is not inspect.Parameter.empty}
    assert defaults["product"] == "3d" and defaults["time_handling"] == "fuzzy"
    assert defaults["interp_type"] == "linear" and defaults["keep_outliers"] is False
