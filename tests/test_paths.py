"""Tests for atmosgw.curate.paths."""

import pytest

from atmosgw.curate import paths
from atmosgw.curate.paths import DataDirNotSetError, local_data_dir


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
