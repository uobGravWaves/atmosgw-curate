"""Instrument loaders: observations in a common format.

``load_limb_sounders``: limb sounders, radio occultation, lidar and sampler
output on a [profile x level] grid (port of the MATLAB ``get_limbsounders``).

``load_nadir_sounders``: nadir-sounder swaths (AIRS, IASI) in the swath format
(port of the MATLAB ``prep_airs_3d``).

``load_time_height``: ground-based [time x level] series at one site (the
ALOMAR lidar) in the time_height format (new code).

``load_trajectories``: in-situ measurements along a path (DWD radiosondes)
in the trajectory format (new code).
"""

from atmosgw.curate.obs.limb import load_limb_sounders
from atmosgw.curate.obs.nadir import load_nadir_sounders
from atmosgw.curate.obs.time_height import load_time_height
from atmosgw.curate.obs.trajectory import load_trajectories

__all__ = ["load_limb_sounders", "load_nadir_sounders", "load_time_height",
           "load_trajectories"]
