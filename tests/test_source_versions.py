"""The data versions the files state: the ``source_versions`` attribute, the
``original_file_version`` coordinate (``file_source``) and the warning when
the files read state different versions (obs/_common.py, the readers).

Synthetic MLS files with the Aura attributes (PGEVersion and DOI in
HDFEOS/ADDITIONAL/FILE_ATTRIBUTES), an ACE file (version in its name) and a
Misc file (none). Real data (skipped if absent): SABER 2020, whose files
state two versions (02.07, 02.08) although every file name says v2.0.
"""

import os
import warnings

import h5py
import numpy as np
import pytest
from scipy.io import savemat

from atmosgw.curate.data import validate
from atmosgw.curate.matlab_time import datenum
from atmosgw.curate.obs import load_limb_sounders
from atmosgw.curate.paths import DataDirNotSetError, dataset_dir


def _write_mls(path, day, version, doi):
    """3 profiles x 5 levels on day (yyyy-mm-dd), stating version and doi."""
    group = "HDFEOS/SWATHS/Temperature"
    with h5py.File(path, "w") as f:
        f[f"{group}/Data Fields/L2gpValue"] = np.full((3, 5), 230.0, dtype=np.float32)
        geo = f"{group}/Geolocation Fields"
        f[f"{geo}/Latitude"] = np.array([-10, 0, 10], dtype=np.float32)
        f[f"{geo}/Longitude"] = np.array([0, 1, 2], dtype=np.float32)
        f[f"{geo}/Pressure"] = np.array([100, 30, 10, 3, 1], dtype=np.float32)
        start = ((np.datetime64(day) - np.datetime64("1993-01-01"))
                 / np.timedelta64(1, "s"))
        f[f"{geo}/Time"] = start + np.array([0.0, 3600, 7200])
        attrs = f.create_group("HDFEOS/ADDITIONAL/FILE_ATTRIBUTES").attrs
        attrs["PGEVersion"] = np.bytes_(version)
        attrs["identifier_product_doi"] = np.bytes_(doi)


@pytest.fixture
def mls(tmp_path):
    (tmp_path / "2026").mkdir()
    _write_mls(tmp_path / "2026" / "MLS-Aura_L2GP-Temperature_v05-04-c01_2026d001.he5",
               "2026-01-01", "v5.04", "10.5067/AURA/MLS/DATA2520")
    _write_mls(tmp_path / "2026" / "MLS-Aura_L2GP-Temperature_v06-04-c01_2026d002.he5",
               "2026-01-02", "v6.04", "10.5067/AURA/MLS/DATA2620")
    return tmp_path


def test_one_version(mls):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ds = load_limb_sounders("2026-01-01", "MLS", file_path=mls, file_source=True)
    assert ds.attrs["source_versions"] == (
        "PGEVersion v5.04 (DOI 10.5067/AURA/MLS/DATA2520) (1 file)")
    assert ds.original_file_version.values.tolist() == [
        "PGEVersion v5.04 (DOI 10.5067/AURA/MLS/DATA2520)"]
    validate(ds, "profile")


def test_mixed_versions_warned_and_counted(mls):
    with pytest.warns(UserWarning, match="2 different data versions"):
        ds = load_limb_sounders(["2026-01-01", "2026-01-02"], "MLS", file_path=mls,
                                file_source=True)
    assert ds.attrs["source_versions"] == (
        "PGEVersion v5.04 (DOI 10.5067/AURA/MLS/DATA2520) (1 file); "
        "PGEVersion v6.04 (DOI 10.5067/AURA/MLS/DATA2620) (1 file)")
    assert ds.original_files.values.tolist() == [
        "MLS-Aura_L2GP-Temperature_v05-04-c01_2026d001.he5",
        "MLS-Aura_L2GP-Temperature_v06-04-c01_2026d002.he5"]
    assert ds.original_file_version.values.tolist()[1].startswith("PGEVersion v6.04")


def _profile_file(path, **extra):
    alt = np.tile(np.arange(10.0, 41.0, 2.0), (2, 1))
    savemat(path, {"Alt": alt, "Lat": np.full_like(alt, 30.0),
                   "Lon": np.full_like(alt, 10.0), "Temp": np.full_like(alt, 230.0),
                   "Time": np.full_like(alt, datenum(2007, 1, 1) + 0.5), **extra})


def test_version_from_file_name_and_none(tmp_path):
    (tmp_path / "2007").mkdir()
    _profile_file(tmp_path / "2007" / "ace_v52_2007d001.mat")
    ds = load_limb_sounders("2007-01-01", "ACE", file_path=tmp_path)
    assert ds.attrs["source_versions"] == "v52 (file name) (1 file)"
    _profile_file(tmp_path / "mine_2007d001.mat")
    ds = load_limb_sounders("2007-01-01", "Misc", misc_path=tmp_path,
                            misc_identifier="mine")
    assert ds.attrs["source_versions"] == "not stated in the files (1 file)"
    ds = load_limb_sounders("2007-01-05", "Misc", misc_path=tmp_path,
                            misc_identifier="mine")
    assert ds.attrs["source_versions"] == ""  # (no file read)


def test_real_saber_2020_states_two_versions():
    try:
        folder = dataset_dir("SABER")
    except (DataDirNotSetError, FileNotFoundError) as error:
        pytest.skip(f"no SABER data: {error}")
    if not os.path.isdir(folder):
        pytest.skip("no SABER data")
    with pytest.warns(UserWarning, match="different data versions"):
        ds = load_limb_sounders(["2020-01-01", "2020-12-31"], "SABER",
                                original_z=True, lat_range=(0, 1))
    assert "Data_Product_Version 02.07" in ds.attrs["source_versions"]
    assert "Data_Product_Version 02.08" in ds.attrs["source_versions"]
