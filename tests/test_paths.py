"""Tests for atmosgw.curate.paths."""

import pytest

from atmosgw.curate import paths
from atmosgw.curate.paths import (
    DataDirNotSetError,
    dataset_dir,
    local_data_dir,
    reference_dir,
)


@pytest.fixture
def clean_env(monkeypatch, tmp_path):
    """No environment variables set, and a home folder with no settings file."""
    monkeypatch.delenv(paths.DATA_DIR_ENV, raising=False)
    monkeypatch.delenv(paths.CONFIG_ENV, raising=False)
    monkeypatch.setattr(paths.Path, "home", lambda: tmp_path / "home")
    return tmp_path


def write_config(home, text):
    config = home / ".config" / "atmosgw" / "config.toml"
    config.parent.mkdir(parents=True)
    config.write_text(text)
    return config


def test_argument(clean_env):
    assert local_data_dir(clean_env) == clean_env


def test_environment_variable(clean_env, monkeypatch):
    monkeypatch.setenv(paths.DATA_DIR_ENV, str(clean_env))
    assert local_data_dir() == clean_env


def test_config_file(clean_env):
    data = clean_env / "data"
    data.mkdir()
    write_config(clean_env / "home", f'data_dir = "{data.as_posix()}"\n')
    assert local_data_dir() == data


def test_config_file_named_by_environment(clean_env, monkeypatch):
    data = clean_env / "data"
    data.mkdir()
    config = clean_env / "elsewhere.toml"
    config.write_text(f'data_dir = "{data.as_posix()}"\n')
    monkeypatch.setenv(paths.CONFIG_ENV, str(config))
    assert local_data_dir() == data


def test_precedence(clean_env, monkeypatch):
    by_arg, by_env, by_file = (clean_env / n for n in ("arg", "env", "file"))
    for folder in (by_arg, by_env, by_file):
        folder.mkdir()
    write_config(clean_env / "home", f'data_dir = "{by_file.as_posix()}"\n')
    assert local_data_dir() == by_file
    monkeypatch.setenv(paths.DATA_DIR_ENV, str(by_env))
    assert local_data_dir() == by_env
    assert local_data_dir(by_arg) == by_arg


def test_tilde_expanded(clean_env, monkeypatch):
    monkeypatch.setenv("HOME", str(clean_env))
    monkeypatch.setenv("USERPROFILE", str(clean_env))  # Windows
    assert local_data_dir("~") == clean_env


def test_nothing_set(clean_env):
    with pytest.raises(DataDirNotSetError, match="data_dir"):
        local_data_dir()


def test_config_without_data_dir(clean_env):
    write_config(clean_env / "home", 'something_else = "x"\n')
    with pytest.raises(DataDirNotSetError):
        local_data_dir()


def test_missing_folder(clean_env):
    missing = clean_env / "not_there"
    with pytest.raises(FileNotFoundError, match="not_there"):
        local_data_dir(missing)
    assert local_data_dir(missing, must_exist=False) == missing


# dataset_dir and reference_dir


@pytest.fixture
def data_root(clean_env):
    """A data root set in the settings file; returns (root, home)."""
    root = clean_env / "data"
    root.mkdir()
    home = clean_env / "home"
    return root, home


def test_dataset_default(data_root):
    root, home = data_root
    write_config(home, f'data_dir = "{root.as_posix()}"\n')
    (root / "SABER" / "raw").mkdir(parents=True)
    assert dataset_dir("SABER") == root / "SABER" / "raw"
    assert dataset_dir("SABER", "pws", must_exist=False) == root / "SABER" / "pws"
    assert dataset_dir("ALOMAR", must_exist=False) == root / "lidar" / "ALOMAR"


def test_dataset_overrides(data_root, clean_env):
    root, home = data_root
    elsewhere = clean_env / "big_disk" / "saber"
    elsewhere.mkdir(parents=True)
    write_config(
        home,
        f'data_dir = "{root.as_posix()}"\n'
        "[datasets.SABER]\n"
        f'raw = "{elsewhere.as_posix()}"\n'
        'pws = "SABER_v2/pws"\n',
    )
    assert dataset_dir("SABER") == elsewhere  # absolute override
    assert dataset_dir("SABER", "pws", must_exist=False) == root / "SABER_v2" / "pws"
    # other datasets keep their defaults
    assert dataset_dir("HIRDLS", must_exist=False) == root / "HIRDLS" / "raw"


def test_dataset_added_by_user(data_root):
    root, home = data_root
    write_config(
        home, f'data_dir = "{root.as_posix()}"\n[datasets.NEWSAT]\nraw = "new/raw"\n'
    )
    assert dataset_dir("NEWSAT", must_exist=False) == root / "new" / "raw"


def test_dataset_unknown(data_root):
    root, home = data_root
    write_config(home, f'data_dir = "{root.as_posix()}"\n')
    with pytest.raises(KeyError, match="Unknown dataset"):
        dataset_dir("saber")  # names are case-sensitive
    with pytest.raises(KeyError, match="no 'pws' folder"):
        dataset_dir("ALOMAR", "pws")


def test_dataset_missing_folder(data_root):
    root, home = data_root
    write_config(home, f'data_dir = "{root.as_posix()}"\n')
    with pytest.raises(FileNotFoundError, match="package defaults"):
        dataset_dir("MLS")


def test_reference_dir(data_root, clean_env):
    root, home = data_root
    (root / "tests").mkdir()
    config = write_config(
        home, f'data_dir = "{root.as_posix()}"\nreference_dir = "tests"\n'
    )
    assert reference_dir() == root / "tests"
    absolute = clean_env / "refs"
    absolute.mkdir()
    config.write_text(f'reference_dir = "{absolute.as_posix()}"\n')
    assert reference_dir() == absolute


def test_reference_dir_unset(data_root):
    root, home = data_root
    write_config(home, f'data_dir = "{root.as_posix()}"\n')
    with pytest.raises(DataDirNotSetError, match="reference_dir"):
        reference_dir()
