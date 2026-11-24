"""Instrument loaders: observations in a common format.

``load_limb_sounders``: limb sounders, radio occultation, lidar and sampler
output on a [profile x level] grid (port of the MATLAB ``get_limbsounders``).
"""

from atmosgw.curate.obs.limb import load_limb_sounders

__all__ = ["load_limb_sounders"]
