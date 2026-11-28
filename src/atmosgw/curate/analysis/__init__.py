"""Analysis steps applied to loader output: detrending and smoothing.

``detrend`` (``analysis.detrend``) removes the background with the group's
current methods; ``smooth_array`` (``analysis.smooth``) is its box or
Gaussian smoothing step; ``find_tropopause_kink`` finds the tropopause
kinks (sharp reversals of the temperature gradient, most likely the
tropopause) that ``detrend`` masks around.
"""

from atmosgw.curate.analysis.detrend import (
    TropopauseKinkWarning,
    detrend,
    find_tropopause_kink,
)
from atmosgw.curate.analysis.smooth import smooth_array

__all__ = ["TropopauseKinkWarning", "detrend", "find_tropopause_kink", "smooth_array"]
