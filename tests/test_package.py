"""Packaging tests: the namespace layout must stay installable alongside other
AtmosGW family members (atmosgw-fuse etc.)."""

import importlib.util

import atmosgw
import atmosgw.curate


def test_version_is_set():
    assert atmosgw.curate.__version__
    assert atmosgw.curate.__version__ != "0.0.0+unknown"


def test_atmosgw_is_a_namespace_package():
    # A namespace package has no __init__.py, so no __file__. If this fails,
    # someone has added src/atmosgw/__init__.py, which would break
    # co-installation with atmosgw-fuse.
    assert getattr(atmosgw, "__file__", None) is None
    assert importlib.util.find_spec("atmosgw").origin in (None, "namespace")
