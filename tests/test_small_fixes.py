"""Small fixes: empty regridded nadir results have the regridded shape;
AIRS 2-D fill values; IASI screen padding for explicit granules follows the
window; the limb ALOMAR reader finds sessions that cross New Year (Q16)."""

import pytest

from atmosgw.curate.obs import load_limb_sounders, load_nadir_sounders, load_time_height
from atmosgw.curate.obs.nadir import _with_neighbours


def test_empty_regridded_swath_has_the_regridded_shape(tmp_path):
    with pytest.warns(UserWarning):
        ds = load_nadir_sounders(["2010-07-01T00:00", "2010-07-01T00:20"], "IASI",
                                 file_path=tmp_path)
    assert ds.attrs["grid"] == "regular"
    assert ds.sizes["cross_track"] == 90 and ds.sizes["along_track"] == 0
    assert "distance_cross" in ds and ds.sizes["level"] == 31


def test_padding_follows_the_window():
    # (the default: half the 10-minute window plus the 1-minute block step)
    assert _with_neighbours([(733043, 1)], 6) == [(733042, 240), (733043, 1),
                                                   (733043, 2)]
    got = _with_neighbours([(733043, 1)], 13)  # (a 26-minute window: 3 each side)
    assert len(got) == 7 and got[0] == (733042, 238) and got[-1] == (733043, 4)


def test_alomar_session_across_new_year(tmp_path):
    from test_time_height_known import sessions_2007, write_year
    write_year(tmp_path, 2007, sessions_2007())
    th = load_time_height("2008-01-01", "ALOMAR", file_path=tmp_path)
    with pytest.warns(UserWarning):  # (the ALOMAR pointer to load_time_height)
        limb = load_limb_sounders("2008-01-01", "ALOMAR", file_path=tmp_path,
                                  original_z=True, file_source=True)
    assert th.sizes["time"] > 0
    assert limb.sizes["profile"] == th.sizes["time"]
    assert limb.original_files.values.tolist() == [th.attrs["files_read"]]
