"""Instrument loaders: observations in a common format.

``load_limb_sounders``: limb sounders, radio occultation, lidar and sampler
output on a [profile x level] grid (port of the MATLAB ``get_limbsounders``).

``load_nadir_sounders``: nadir-sounder swaths (AIRS) in the swath format
(port of the MATLAB ``prep_airs_3d``).
"""

from atmosgw.curate.obs.limb import load_limb_sounders
from atmosgw.curate.obs.nadir import load_nadir_sounders

__all__ = ["load_limb_sounders", "load_nadir_sounders"]
