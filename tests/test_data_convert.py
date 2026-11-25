"""Tests for atmosgw.curate.data.convert.to_profiles (synthetic AIRS swath).

Exact checks: shapes, order, values (reshaped, not computed) and states; the
estimated pressure is checked against atmosphere.height_to_pressure itself.
"""

import numpy as np
import pytest
from test_nadir_known import D, seconds, write_granule

from atmosgw.curate.atmosphere import height_to_pressure
from atmosgw.curate.data import states, to_profiles, validate
from atmosgw.curate.obs import load_nadir_sounders


@pytest.fixture
def swath(tmp_path):
    temp = np.broadcast_to(220 + 0.5 * np.r_[0:61:3, 65:91:5],
                           (135, 90, 27)).astype(np.float32).copy()
    temp[100, 20, 4] = -999.0  # missing
    write_granule(tmp_path, 53, seconds(53), temp)
    return load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=tmp_path,
                               n_cross=12, n_along=9, day_night=True)


def test_shapes_order_and_values(swath):
    p = to_profiles(swath)
    validate(p, "profile")
    n_al, n_x, n_lev = (swath.sizes[d] for d in ("along_track", "cross_track", "level"))
    assert p.sizes["profile"] == n_al * n_x and p.sizes["level"] == n_lev
    k = 4 * n_x + 7  # along 5, cross 8 (1-based)
    assert p.swath_along.values[k] == 5 and p.swath_cross.values[k] == 8
    np.testing.assert_array_equal(p.temperature.values[k],
                                  swath.temperature.values[4, 7])
    np.testing.assert_array_equal(p.temperature_flag.values[k],
                                  swath.temperature_flag.values[4, 7])
    assert np.all(p.latitude.values[k] == swath.latitude.values[4, 7])
    np.testing.assert_array_equal(p.altitude.values[k], swath.altitude.values)
    assert p.granule.values[k] == 53
    assert p.day_retrieval.values[k] == swath.day_retrieval.values[4, 7]
    assert p.attrs["converted_from"] == "swath" and p.attrs["format"] == "profile"


def test_pressure_is_derived(swath):
    p = to_profiles(swath)
    np.testing.assert_array_equal(p.pressure.values,
                                  np.asarray(height_to_pressure(p.altitude.values),
                                             dtype=np.float32))
    assert p.pressure.attrs["origin"] == "derived"
    assert states.has_state(p.pressure_flag.values, "derived").all()


def test_file_provenance_round_trip(tmp_path):
    for g in (53, 54):
        write_granule(tmp_path, g, seconds(g))
    sw = load_nadir_sounders(None, "AIRS", granules=[(D, 53), (D, 54)],
                             file_path=tmp_path, n_cross=12, n_along=9,
                             file_source=True)
    p = to_profiles(sw)
    validate(p, "profile")
    n_x, n_lev = sw.sizes["cross_track"], sw.sizes["level"]
    assert p.source_file.dims == ("profile", "level")
    assert p.source_row.dims == ("profile",)
    assert p.source_file.dtype == sw.source_file.dtype
    assert list(p.original_files.values) == list(sw.original_files.values)
    assert sorted(set(sw.source_file.values)) == [1, 2]  # both files used
    for row in (0, 4, sw.sizes["along_track"] - 1):
        for cross in (0, 7, n_x - 1):
            k = row * n_x + cross
            assert p.swath_along.values[k] == row + 1
            assert np.all(p.source_file.values[k] == sw.source_file.values[row])
            assert p.source_row.values[k] == sw.source_row.values[row]
    np.testing.assert_array_equal(
        p.source_file.values,
        np.repeat(np.repeat(sw.source_file.values, n_x)[:, None], n_lev, axis=1))
    avail = states.parse_availability(p.attrs["availability"])
    assert avail["source_profile"] == "not_provided"
    assert "source_file" not in avail and "original_files" not in avail


def test_availability_inherited(swath):
    avail = states.parse_availability(to_profiles(swath).attrs["availability"])
    assert avail["source_file"] == "not_requested"
    assert avail["original_files"] == "not_requested"
    assert avail["source_profile"] == "not_provided"
    text = swath.attrs["availability"].replace("source_file: not_requested",
                                               "source_file: not_in_files")
    p = to_profiles(swath.assign_attrs(availability=text))
    validate(p, "profile")
    avail = states.parse_availability(p.attrs["availability"])
    assert avail["source_file"] == "not_in_files"


def test_rejects_2d_and_other_formats(swath):
    with pytest.raises(ValueError):
        to_profiles(swath.drop_vars(["temperature", "temperature_flag", "altitude"]))
    with pytest.raises(ValueError):
        to_profiles(swath.assign_attrs(format="profile"))
