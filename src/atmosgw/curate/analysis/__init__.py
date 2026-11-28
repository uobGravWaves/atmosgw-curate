"""Analysis steps applied to loader output: detrending and smoothing.

``detrend`` (``analysis.detrend``) removes the background with the group's
current methods; ``smooth_array`` (``analysis.smooth``) is its box or
Gaussian smoothing step; ``find_tropopause`` gives the lapse-rate
tropopause that ``detrend`` masks around.
"""

from atmosgw.curate.analysis.detrend import (
    TropopauseMaskWarning,
    detrend,
    find_tropopause,
)
from atmosgw.curate.analysis.smooth import smooth_array

__all__ = ["TropopauseMaskWarning", "detrend", "find_tropopause", "smooth_array"]
