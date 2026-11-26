"""
Files need not be in the usual year (and day) folders: each reader looks for
a file by name anywhere below the dataset folder if it is not at its usual
path. Each test writes the same files in the usual layout and flat (or in
other folders), and checks the loads are the same.
"""

import shutil
import warnings

import numpy as np
import pytest
import xarray as xr
from scipy.io import savemat
from test_limb_gnss_ace import D1, _ace, _gnss_content
from test_nadir_airs2d import DAY as DAY_2D
from test_nadir_airs2d import write as write_2d
from test_nadir_iasi import DAY_SECONDS, write_orbit
from test_nadir_known import D, seconds, write_granule

from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs import _limb_ace_gnss, load_limb_sounders, load_nadir_sounders


def flatten(usual, flat, nest=""):
    """Copy every file below usual into flat/nest (no year or day folders)."""
    target = flat / nest
    target.mkdir(parents=True, exist_ok=True)
    for path in usual.rglob("*"):
        if path.is_file():
            shutil.copy(path, target / path.name)
    return flat


def same(a, b):
    """The same data and coordinates (attributes may name other folders)."""
    xr.testing.assert_equal(a, b)
    assert a.attrs["files_read"] == b.attrs["files_read"]
    assert a.attrs.get("files_not_found", "") == b.attrs.get("files_not_found", "")


@pytest.mark.parametrize("nest", ["", "some/other/folder"])
def test_ace_flat(tmp_path, nest):
    (tmp_path / "usual").mkdir()
    _ace(tmp_path / "usual")
    flat = flatten(tmp_path / "usual", tmp_path / "flat", nest)
    same(load_limb_sounders(D1, "ACE", file_path=tmp_path / "usual"),
         load_limb_sounders(D1, "ACE", file_path=flat))


def test_gnss_flat(tmp_path, monkeypatch):
    (tmp_path / "usual" / "2007").mkdir(parents=True)
    (tmp_path / "usual" / "2007" / "merged_ro_2007_001.mat").write_bytes(b"")
    flat = flatten(tmp_path / "usual", tmp_path / "flat")
    read = []
    monkeypatch.setattr(_limb_ace_gnss, "read_mat",
                        lambda path: read.append(path) or _gnss_content())
    same(load_limb_sounders(D1, "GNSS", file_path=tmp_path / "usual"),
         load_limb_sounders(D1, "GNSS", file_path=flat))
    assert read[1] == str(flat / "merged_ro_2007_001.mat")


def test_pw_flat(tmp_path):
    lev = np.arange(10.0, 61.0, 5.0)
    alt = 0.3 * np.arange(1, 5)[:, None] + lev
    content = {"Alt": alt, "Lat": np.full_like(alt, 20.0),
               "Lon": np.full_like(alt, 5.0),
               "Time": np.full_like(alt, datenum(2007, 1, 1) + 0.5),
               "Temp_Residual": 0.1 * alt,
               "Temp_PW": alt[:, :, None] * np.arange(1, 4)}
    (tmp_path / "usual" / "2007").mkdir(parents=True)
    savemat(tmp_path / "usual" / "2007" / "SABER_PWs_2007d001.mat", content)
    flat = flatten(tmp_path / "usual", tmp_path / "flat")
    kw = {"hindley_pws": True}
    same(load_limb_sounders(D1, "SABER", file_path=tmp_path / "usual", **kw),
         load_limb_sounders(D1, "SABER", file_path=flat, **kw))


def test_airs_flat_and_duplicates(tmp_path):
    for g in (53, 54):
        write_granule(tmp_path / "usual", g, seconds(g))
    flat = flatten(tmp_path / "usual", tmp_path / "flat")
    usual = load_nadir_sounders(None, "AIRS", granules=[(D, 53), (D, 54)],
                                file_path=tmp_path / "usual")
    same(usual, load_nadir_sounders(None, "AIRS", granules=[(D, 53), (D, 54)],
                                    file_path=flat))
    # the same file name in two folders: the first (search order) is used, and
    # both are named, by their paths below the dataset folder, and warned
    flatten(tmp_path / "usual", flat, "copy")
    with pytest.warns(UserWarning, match="more than one matching"):
        twice = load_nadir_sounders(None, "AIRS", granules=[(D, 53)], file_path=flat)
    assert "copy/airs_2007_001_053.nc" in twice.attrs["multiple_matches"]


def test_airs_usual_layout_wins(tmp_path):
    # a file at its usual path is used even if a file of that name is elsewhere
    write_granule(tmp_path / "usual", 53, seconds(53))
    flatten(tmp_path / "usual", tmp_path / "usual", "elsewhere")
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # (no multiple-match warning)
        ds = load_nadir_sounders(None, "AIRS", granules=[(D, 53)],
                                 file_path=tmp_path / "usual")
    assert ds.attrs.get("multiple_matches", "") == ""


def test_airs2d_flat(tmp_path):
    write_2d(tmp_path / "usual", 121, year_folder=True)
    flat = flatten(tmp_path / "usual", tmp_path / "flat")
    kw = {"product": "2d", "granules": [(DAY_2D, 121)]}
    same(load_nadir_sounders(None, "AIRS", file_path=tmp_path / "usual", **kw),
         load_nadir_sounders(None, "AIRS", file_path=flat, **kw))


def test_iasi_flat_uses_the_date_in_the_name(tmp_path):
    write_orbit(tmp_path / "usual", DAY_SECONDS, int(3056 / 4))
    write_orbit(tmp_path / "usual", DAY_SECONDS + 3056, 900)
    flat = flatten(tmp_path / "usual", tmp_path / "flat")
    write_orbit(tmp_path / "other", DAY_SECONDS + 86400 * 400, 900)  # another date
    flatten(tmp_path / "other", flat)
    kw = {"granules": [("2007-11-07", 9)], "original_grid": True}
    a = load_nadir_sounders(None, "IASI", file_path=tmp_path / "usual", **kw)
    b = load_nadir_sounders(None, "IASI", file_path=flat, **kw)
    same(a, b)
    assert b.attrs["files_read"].count("iasi_") == 2
    # a day with no files of its own, in a flat folder, is recorded as missing
    with pytest.warns(UserWarning, match="2007-11-08"):
        t = load_nadir_sounders(("2007-11-08", "2007-11-08 00:20"), "IASI",
                                file_path=flat, original_grid=True)
    assert "2007-11-08" in t.attrs["files_not_found"]
