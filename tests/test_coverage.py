"""Loads that return less than was asked for say so: a missing SABER month
in the middle of a month; files read that hold nothing in the requested
times (limb, nadir, time-height); an IASI day split between folder layouts,
and a day folder with no file of that day. Nadir scan widths are judged
against the instrument's nominal width, so a collapsed granule is found
alone and does not affect a healthy one loaded beside it (Q14). Synthetic
files only, so these run in CI."""

import shutil
import warnings

import numpy as np
import pytest
from scipy.io import savemat
from test_nadir_geometry import _granule
from test_nadir_iasi import DAY_SECONDS, write_orbit
from test_nadir_known import D, seconds
from test_readers_synthetic import _saber
from test_time_height_known import sessions_2007, write_year

from atmosgw.curate.data import states, validate
from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs import load_limb_sounders, load_nadir_sounders, load_time_height


def _messages(call):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = call()
    return out, [str(w.message) for w in caught]


# -------------------------------------------------- nothing in the times asked for


def test_saber_month_missing_mid_month_is_reported(tmp_path):
    _saber(tmp_path)  # (January 2007 only)
    ds, said = _messages(lambda: load_limb_sounders("2007-02-15", "SABER",
                                                    file_path=tmp_path))
    assert ds.sizes["profile"] == 0
    assert ds.attrs["files_not_found"] == "February2007"
    assert any("February2007" in m for m in said)
    # a range starting mid-month: the month is reported once
    ds, said = _messages(lambda: load_limb_sounders(["2007-02-15", "2007-02-20"],
                                                    "SABER", file_path=tmp_path))
    assert ds.attrs["files_not_found"] == "February2007"
    # the look-back day alone is not reported (1 February: January is read)
    ds, _ = _messages(lambda: load_limb_sounders("2007-02-01", "SABER",
                                                 file_path=tmp_path))
    assert ds.attrs["files_not_found"] == "February2007"


def _misc_file(folder, day_of_data):
    alt = np.tile(np.arange(10.0, 41.0, 2.0), (2, 1))
    savemat(folder / "mine_2007d001.mat", {
        "Alt": alt, "Lat": np.full_like(alt, 30.0), "Lon": np.full_like(alt, 10.0),
        "Temp": np.full_like(alt, 230.0),
        "Time": np.full_like(alt, datenum(2007, 1, day_of_data) + 0.5)})
    return dict(misc_path=folder, misc_identifier="mine")


def test_limb_file_holding_nothing_in_the_times_warns(tmp_path):
    opts = _misc_file(tmp_path, 5)  # (named 1 January, holds 5 January)
    ds, said = _messages(lambda: load_limb_sounders("2007-01-01", "Misc", **opts))
    assert ds.sizes["profile"] == 0
    found = [m for m in said if "hold no profiles in the requested times" in m]
    assert len(found) == 1 and "2007-01-05T12:00:00" in found[0]
    validate(ds, "profile")
    # silenced by date_warning=False; no such warning when there are data
    _, said = _messages(lambda: load_limb_sounders("2007-01-01", "Misc",
                                                   date_warning=False, **opts))
    assert not said
    opts = _misc_file(tmp_path, 1)
    ds, said = _messages(lambda: load_limb_sounders("2007-01-01", "Misc", **opts))
    assert ds.sizes["profile"] == 2 and not any("hold no" in m for m in said)


def test_limb_box_with_nothing_inside_does_not_warn(tmp_path):
    opts = _misc_file(tmp_path, 1)
    ds, said = _messages(lambda: load_limb_sounders("2007-01-01", "Misc",
                                                    lat_range=(-10, 0), **opts))
    assert ds.sizes["profile"] == 0 and not any("hold no" in m for m in said)


def test_nadir_request_between_granule_mean_times_warns(tmp_path):
    _granule(tmp_path, 53)
    start = np.datetime64("2000-01-01") + np.timedelta64(int(seconds(53)), "s")
    # (the granule's mean time is about 3 minutes after its start)
    times = [str(start + np.timedelta64(30, "s")), str(start + np.timedelta64(90, "s"))]
    ds, said = _messages(lambda: load_nadir_sounders(times, "AIRS", file_path=tmp_path,
                                                     time_handling="strict"))
    assert ds.sizes["along_track"] == 0
    found = [m for m in said if "hold no granules" in m]
    assert len(found) == 1 and "The 1 AIRS file(s) read" in found[0]
    assert not any("No AIRS file was read" in m for m in said)
    validate(ds, "swath")


def test_nadir_box_with_nothing_inside_does_not_warn(tmp_path):
    _granule(tmp_path, 53)
    ds, said = _messages(lambda: load_nadir_sounders(
        None, "AIRS", file_path=tmp_path, granules=[(D, 53)], lat_range=(-80, -70)))
    assert ds.sizes["along_track"] == 0 and not any("hold no" in m for m in said)


def test_single_strict_instant_warns_once(tmp_path):
    opts = _misc_file(tmp_path, 1)
    ds, said = _messages(lambda: load_limb_sounders(
        "2007-01-01T03:00", "Misc", time_handling="strict", **opts))
    assert ds.sizes["profile"] == 0
    assert any("exactly that instant" in m for m in said)
    assert not any("hold no" in m for m in said)


def test_time_height_day_without_a_session_warns(tmp_path):
    write_year(tmp_path, 2007, sessions_2007())
    ds, said = _messages(lambda: load_time_height("2007-01-01", "ALOMAR",
                                                  file_path=tmp_path))
    assert ds.sizes["time"] == 0
    assert any("hold no measurements in the requested times" in m for m in said)
    ds, said = _messages(lambda: load_time_height("2007-01-23", "ALOMAR",
                                                  file_path=tmp_path))
    assert ds.sizes["time"] > 0 and not any("hold no" in m for m in said)


# -------------------------------------------------- IASI folder layouts


def _iasi_day(root):
    """Two orbit files of 2007-11-07, 00:00 and 01:00, in 2007/311."""
    a = write_orbit(root, DAY_SECONDS + 60, 90)
    b = write_orbit(root, DAY_SECONDS + 3600, 90)
    return a, b


def test_iasi_day_split_between_layouts_is_read_whole(tmp_path):
    a, b = _iasi_day(tmp_path)
    whole = load_nadir_sounders("2007-11-07", "IASI", file_path=tmp_path,
                                screen="none", original_grid=True)
    # the second orbit moved out of the day folder, to the top of the tree
    shutil.move(b, tmp_path / b.name)
    split, said = _messages(lambda: load_nadir_sounders(
        "2007-11-07", "IASI", file_path=tmp_path, screen="none", original_grid=True))
    assert split.sizes["along_track"] == whole.sizes["along_track"] == 180
    np.testing.assert_array_equal(split.temperature.values, whole.temperature.values)
    # the day folder empty, both files elsewhere
    shutil.move(a, tmp_path / a.name)
    moved, _ = _messages(lambda: load_nadir_sounders(
        "2007-11-07", "IASI", file_path=tmp_path, screen="none", original_grid=True))
    assert moved.sizes["along_track"] == 180


def test_iasi_day_folder_without_its_files_is_reported(tmp_path):
    _iasi_day(tmp_path)
    (tmp_path / "2007" / "312").mkdir()  # (8 November: a folder, no file)
    ds, said = _messages(lambda: load_nadir_sounders(
        ["2007-11-07", "2007-11-08"], "IASI", file_path=tmp_path, screen="none",
        original_grid=True))
    assert "2007-11-08" in ds.attrs["files_not_found"]
    assert any("not found" in m for m in said)


def test_iasi_misfiled_file_in_another_days_folder_is_left_out(tmp_path):
    a, _ = _iasi_day(tmp_path)
    wrong = tmp_path / "2007" / "149"
    wrong.mkdir(parents=True)
    other = write_orbit(tmp_path, DAY_SECONDS + 7200, 90)  # (02:00 on 7 November)
    shutil.move(other, wrong / other.name)  # (filed under 29 May)
    ds, _ = _messages(lambda: load_nadir_sounders(
        "2007-11-07", "IASI", file_path=tmp_path, screen="none", original_grid=True))
    assert ds.sizes["along_track"] == 180
    assert other.name not in ds.attrs["files_read"]


def test_iasi_copies_on_a_day_not_used_are_not_reported(tmp_path):
    _iasi_day(tmp_path)
    nxt = write_orbit(tmp_path, DAY_SECONDS + 86400 + 7200, 90)  # (8 Nov, 02:00)
    (tmp_path / "copies").mkdir()
    shutil.copy(nxt, tmp_path / "copies" / nxt.name)
    ds, said = _messages(lambda: load_nadir_sounders(
        ["2007-11-07T00:00", "2007-11-07T02:00"], "IASI", file_path=tmp_path,
        screen="none", original_grid=True))
    assert ds.attrs["multiple_matches"] == ""
    assert not any("more than one matching" in m for m in said)


def test_iasi_same_file_in_two_places_is_read_once(tmp_path):
    a, _ = _iasi_day(tmp_path)
    (tmp_path / "copies").mkdir()
    shutil.copy(a, tmp_path / "copies" / a.name)
    ds, said = _messages(lambda: load_nadir_sounders(
        "2007-11-07", "IASI", file_path=tmp_path, screen="none", original_grid=True))
    assert ds.sizes["along_track"] == 180
    assert a.name in ds.attrs["multiple_matches"]
    assert any("1 day(s) had more than one matching" in m for m in said)


# -------------------------------------------------- scan widths (Q14)


@pytest.mark.parametrize("order", [(53, 54), (54, 53)])
def test_collapsed_granule_alone_or_beside_a_healthy_one(tmp_path, order):
    _granule(tmp_path, 53)
    _granule(tmp_path, 54, cross_scale=0.01)  # (scans 1 % of the usual width)

    def load(granules):
        return _messages(lambda: load_nadir_sounders(
            None, "AIRS", file_path=tmp_path, granules=[(D, g) for g in granules]))

    healthy, _ = load([53])
    assert healthy.attrs["grid"] == "regular"
    assert np.isfinite(healthy.temperature.values).all()
    collapsed, said = load([54])  # (found alone: it cannot be placed)
    assert collapsed.attrs["grid"] == "native"
    assert np.all(states.has_state(collapsed.temperature_flag.values,
                                   "bad_geolocation"))
    assert any("scan widths" in m for m in said)
    both, _ = load(order)
    assert both.attrs["grid"] == "regular"
    assert set(np.unique(both.granule.values)) == {53}  # (54 left out, as missing)
    for name in ("temperature", "latitude", "longitude", "distance_cross"):
        np.testing.assert_array_equal(both[name].values, healthy[name].values)
    validate(both, "swath")
