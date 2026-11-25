"""Common data formats, data states, and their checker.

``validate(ds, format_name)``: check a Dataset against a format (``"swath"``,
``"profile"``); ``availability(ds, name)``: whether a variable is present, or
why not. See ``formats.py``, ``states.py`` and
``docs/design/data_formats.md``.
"""

from atmosgw.curate.data.formats import (
    FORMATS,
    FormatError,
    availability,
    problems,
    validate,
)

__all__ = ["FORMATS", "FormatError", "availability", "problems", "validate"]
