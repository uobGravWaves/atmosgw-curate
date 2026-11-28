"""
Data states: why a value is absent, or what kind of value it is.

The shared registry of the states recorded in every loader's output, and
helpers to set, test and describe them. Two levels (see
``docs/design/data_formats.md`` section 6):

- **Whole variables** (``availability``): a variable the format defines but
  the Dataset does not hold is listed in the dataset attribute
  ``availability`` with one of ``AVAILABILITY`` (e.g. ``not_in_files``), so
  "never available" is never confused with "missing" or "bad".
- **Individual points** (``STATES``): each data variable ``<name>`` has an
  unsigned-integer variable ``<name>_flag`` of the same shape, one bit per
  state, following the CF flag convention (``flag_masks``,
  ``flag_meanings``). Zero: none of the recorded states applies.

Inputs and outputs
------------------
================================================  ==============================
Function                                          Returns
================================================  ==============================
``state(name)``                                   the ``State`` called ``name``
``flag_attrs()``                                  CF attributes for a flag
                                                  variable (dict)
``new_flags(shape)``                              zeros of ``FLAG_DTYPE``
``set_state(flags, where, name)``                 ``flags`` with the bit set
                                                  where ``where`` is True
                                                  (in place; also returned)
``has_state(flags, name)``                        bool array
``no_value_mask(flags)``                          bool array: any state whose
                                                  value is NaN
``state_counts(flags)``                           {state name: count}, non-zero
                                                  only
``format_availability(mapping)``                  attribute text
``parse_availability(text)``                      {variable: state}
================================================  ==============================

New states are added by appending to ``STATES`` with the next free bit;
nothing else changes, and saved files stay readable because each carries
its own ``flag_masks`` and ``flag_meanings``.

The thresholds behind some states (e.g. 100-400 K) belong to the code that
sets them, not to this registry.

Origin
------
New code (no MATLAB equivalent): the MATLAB loaders set bad points to NaN
without recording why.

Validation
----------
``tests/test_data_formats.py`` (the registry: unique bits, setting,
testing and counting states, CF attributes, availability text);
``tests/test_states_conformance.py`` (the states each loader sets).

Status
------
Proposal (the state list is a draft for agreement). AI-originated
(Claude), from the design agreed with C. Wright on 2026-11-25; to be
reviewed by C. Wright.

Usage
-----
>>> import numpy as np
>>> from atmosgw.curate.data.states import new_flags, set_state, has_state
>>> f = new_flags(3)
>>> _ = set_state(f, np.array([True, False, False]), "missing")
>>> has_state(f, "missing").tolist()
[True, False, False]

Change log
----------
2026/11/25 [Claude] initial version
2026/11/25 [Claude] descriptions match the shared rules; counts_text
2026/11/26 [Claude] docstring: Validation names the existing tests
2026/11/27 [Claude] descriptions: trajectories; instrument_flag rules referenced
2026/11/28 [Claude] qc_incomplete (bit 2048), agreed with C. Wright

"""

from collections import namedtuple

import numpy as np

State = namedtuple("State", "name bit no_value group description")

# the registry: bit values are fixed once released; append new states only
STATES = (
    State("no_observation", 1, True, "never there",
          "nothing measured here: outside the swath, above or below a "
          "profile's measured values, padding between profiles, or a gap "
          "between segments"),
    State("missing", 2, True, "never there",
          "a gap inside what was measured: a fill value or NaN at a scanned "
          "point, or between a profile's first and last value"),
    State("outside_requested_range", 4, True, "never there",
          "removed because outside the height, latitude or longitude range "
          "asked for"),
    State("bad_geolocation", 8, True, "removed as bad",
          "latitude or longitude out of range (nadir: also NaN, time missing, or "
          "a scan of implausible width; limb and trajectories: out of range "
          "only)"),
    State("outside_physical_range", 16, True, "removed as bad",
          "temperature < 100 K (IASI: <= 100 K, its retrieval's lower bound), "
          "or > 400 K (limb, time-height and trajectories: only below 100 km; "
          "nadir: at any height)"),
    State("statistical_outlier", 32, True, "removed as bad",
          "marked by the instrument's statistical screen (obs._screen)"),
    State("instrument_flag", 64, True, "removed as bad",
          "the instrument's own quality control marks it bad (which rules: "
          "docs/design/data_formats.md 6.2)"),
    State("derived", 128, False, "qualified value",
          "computed, not measured or retrieved (e.g. pressure from height)"),
    State("filled", 256, False, "qualified value",
          "filled from neighbours after loading"),
    State("extrapolated", 512, False, "qualified value",
          "outside the source data's coverage (only on values that exist)"),
    State("a_priori_dominated", 1024, False, "qualified value",
          "the retrieval has little information here, as the provider marks "
          "it (HIRDLS: negative precision)"),
    State("qc_incomplete", 2048, False, "qualified value",
          "the provider's quality control could not be completed for this "
          "value: a quality field it needs is missing from the file (MLS); "
          "the other rules were applied"),
)

FLAG_DTYPE = np.uint16  # 16 states; change to uint32 if more are ever needed

# whole-variable states, for the dataset attribute "availability"
AVAILABILITY = {
    "not_measured": "the instrument cannot provide it",
    "not_in_files": "it exists for this instrument, but not in the files read",
    "not_requested": "not asked for in this call",
    "not_provided": "the data supply no per-point value",
}

_BY_NAME = {s.name: s for s in STATES}


def state(name):
    """The State called name (KeyError with the known names if unknown)."""
    try:
        return _BY_NAME[name]
    except KeyError:
        known = ", ".join(_BY_NAME)
        raise KeyError(f"Unknown state {name!r}; known: {known}") from None


def flag_attrs():
    """CF attributes for a flag variable."""
    return {
        "flag_masks": np.array([s.bit for s in STATES], dtype=FLAG_DTYPE),
        "flag_meanings": " ".join(s.name for s in STATES),
        "long_name": "data state (see atmosgw.curate.data.states)",
    }


def new_flags(shape):
    """A flag array of zeros (no state set)."""
    return np.zeros(shape, dtype=FLAG_DTYPE)


def set_state(flags, where, name):
    """Set the bit for state name where `where` is True (in place); returns flags."""
    flags[np.asarray(where, dtype=bool)] |= FLAG_DTYPE(state(name).bit)
    return flags


def has_state(flags, name):
    """True where the bit for state name is set."""
    return (np.asarray(flags) & FLAG_DTYPE(state(name).bit)) != 0


def no_value_mask(flags):
    """True where any state whose value is NaN is set."""
    bits = FLAG_DTYPE(sum(s.bit for s in STATES if s.no_value))
    return (np.asarray(flags) & bits) != 0


def state_counts(flags):
    """{state name: number of points with it}, for states that occur."""
    flags = np.asarray(flags)
    counts = {s.name: int(np.count_nonzero(flags & FLAG_DTYPE(s.bit)))
              for s in STATES}
    return {k: v for k, v in counts.items() if v}


def counts_text(flags_by_variable):
    """The ``state_counts`` attribute text for {variable: flag array}, e.g.
    'temperature: missing 3, derived 10; pressure: none'."""
    parts = []
    for name, flags in flags_by_variable.items():
        found = state_counts(flags)
        parts.append(f"{name}: " + (", ".join(f"{k} {n}" for k, n in found.items())
                                    or "none"))
    return "; ".join(parts)


def format_availability(mapping):
    """Attribute text for {variable: availability state}, e.g. 'a: not_in_files'."""
    for name, value in mapping.items():
        if value not in AVAILABILITY:
            raise ValueError(f"{name}: unknown availability {value!r}; "
                             f"known: {', '.join(AVAILABILITY)}")
    return "; ".join(f"{k}: {v}" for k, v in mapping.items())


def parse_availability(text):
    """{variable: availability state} from the attribute text ('' gives {})."""
    out = {}
    for part in filter(None, (p.strip() for p in str(text).split(";"))):
        name, _, value = (x.strip() for x in part.partition(":"))
        out[name] = value
    return out
