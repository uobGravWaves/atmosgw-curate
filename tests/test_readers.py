"""Tests for atmosgw.curate.readers (find_files, read_mat).

find_files reproduces the MATLAB wildcardsearch/regexpdir rules (see the module
docstring); expected results below follow those rules.
"""

import os

import numpy as np
from scipy.io import savemat

from atmosgw.curate.readers import find_files, read_hdf5, read_mat, read_netcdf


def make_tree(root):
    """root/a.mat, root/B_2007d001.mat, root/sub/x_2007d001.mat, root/sub/deeper/y.txt,
    and a folder root/dir_2007d001/ holding z.mat."""
    for path in ("a.mat", "B_2007d001.mat", "sub/x_2007d001.mat", "sub/deeper/y.txt",
                 "dir_2007d001/z.mat"):
        file = root / path
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_text("")


def rel(paths, root):
    return [os.path.relpath(p, root) + ("/" if p.endswith(os.sep) else "")
            for p in paths]


def test_order_folders_first_then_files(tmp_path):
    make_tree(tmp_path)
    got = rel(find_files(tmp_path, "*.mat"), tmp_path)
    # sub-folders (each searched in turn) before the files at each level
    assert got == ["dir_2007d001/z.mat", "sub/x_2007d001.mat", "B_2007d001.mat",
                   "a.mat"]


def test_matches_relative_path_and_folders(tmp_path):
    make_tree(tmp_path)
    got = rel(find_files(tmp_path, "*2007d001*"), tmp_path)
    # the matching folder itself is listed, as in MATLAB
    assert got == ["dir_2007d001/", "dir_2007d001/z.mat", "sub/x_2007d001.mat",
                   "B_2007d001.mat"]


def test_case_and_wildcards(tmp_path):
    make_tree(tmp_path)
    assert rel(find_files(tmp_path, "b_2007*"), tmp_path) == ["B_2007d001.mat"]
    assert find_files(tmp_path, "b_2007*", case_sensitive=True) == []
    # '?' is one character; the search is not anchored, as in MATLAB
    assert rel(find_files(tmp_path, "x?2007"), tmp_path) == ["sub/x_2007d001.mat"]
    assert len(find_files(tmp_path, "?.mat")) == 4
    # '.' is a literal dot
    assert find_files(tmp_path, "a_mat") == []


def test_alternatives(tmp_path):
    make_tree(tmp_path)
    got = rel(find_files(tmp_path, f"a.mat{os.pathsep}y.txt"), tmp_path)
    assert got == ["sub/deeper/y.txt", "a.mat"]


def test_missing_root(tmp_path):
    assert find_files(tmp_path / "nothing", "*") == []


def test_read_netcdf_as_nph_getnet(tmp_path):
    import netCDF4

    file = tmp_path / "t.nc"
    with netCDF4.Dataset(file, "w") as f:
        f.createDimension("event", 3)
        f.createDimension("level", 2)
        v = f.createVariable("packed", "i2", ("event", "level"), fill_value=-1)
        v.missing_value = np.int16(-999)
        v.scale_factor = 0.5
        v.add_offset = 10.0
        v.valid_range = np.array([0, 5], dtype="i2")  # ignored, as in nph_getnet
        f.set_auto_maskandscale(False)
        v[:] = np.array([[2, -1], [-999, 8], [4, 6]], dtype="i2")
        w = f.createVariable("plain", "f4", ("event",))
        w[:] = np.array([1.5, 2.5, 3.5], dtype="f4")
    d = read_netcdf(file, ["packed", "plain"])
    np.testing.assert_array_equal(d["packed"], [[11.0, np.nan], [np.nan, 14.0],
                                                [12.0, 13.0]])
    assert d["packed"].dtype == np.float64 and d["plain"].dtype == np.float64
    np.testing.assert_array_equal(read_netcdf(file, ["plain"], rows=np.array([0, 2]))
                                  ["plain"], [1.5, 3.5])


def test_read_hdf5_as_double(tmp_path):
    import h5py

    file = tmp_path / "t.h5"
    with h5py.File(file, "w") as f:
        f["g/x"] = np.arange(6, dtype="f4").reshape(2, 3)
    x = read_hdf5(file, "g/x")
    assert x.dtype == np.float64 and x.shape == (2, 3)


def test_read_mat_shapes_and_structs(tmp_path):
    file = tmp_path / "t.mat"
    savemat(file, {"row": np.arange(3.0)[None, :], "col": np.arange(3.0)[:, None],
                   "scalar": 7.0, "cube": np.arange(24.0).reshape(2, 3, 4),
                   "S": {"a": np.ones((2, 2)), "b": 5.0}})
    d = read_mat(file)
    assert d["row"].shape == (1, 3) and d["col"].shape == (3, 1)
    assert d["scalar"].shape == (1, 1)
    np.testing.assert_array_equal(d["cube"], np.arange(24.0).reshape(2, 3, 4))
    assert isinstance(d["S"], dict) and d["S"]["a"].shape == (2, 2)
    assert not any(k.startswith("__") for k in d)
