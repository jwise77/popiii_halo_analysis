################################################################################
#
#    Saves arrays of redshifts and cosmological times corresponding to each
#    dataset for later use.
#
################################################################################

import os
import yt
import ytree
import numpy as np
from file_locations import *

# Load datasets and arbor
ts = yt.load(param_files)
a = ytree.load(tree_file)

# Determine number of datasets.
N = len(ts)

# Loop through halos in the arbor to find the progenitor line with the longest history
# This will be the root halo whose progenitor line will be used as reference for redshifts.
maxN = -1
ref_halo = None
for root_halo in a:
    if len(list(root_halo['prog'])) > maxN:
        maxN = len(list(root_halo['prog']))
        ref_halo = root_halo
print("Storing %d of %d redshifts that have halo data" % (maxN, N))
N0 = N
N = min(N, maxN)

# There are slight discrepancies between redshifts stored in the datasets and
# redshifts in the saved ytree arbor so save those as separate arrays.
ytree_redshifts = np.zeros(N)
ds_redshifts = np.zeros(N)
cosmological_times = np.zeros(N)

# Loop through each dataset and halo in the progenitor line to save redshifts
# and cosmological time.
# Reverse order of progenitor line to go in increasing time.
start = N0-N
for i, ds in enumerate(ts[start:]):
    ytree_redshifts[-(i+1)] = ref_halo['prog', 'redshift'][i]
    ds_redshifts[i] = ds.current_redshift
    cosmological_times[i] = ds.current_time.to('yr')

# Save all arrays to directory for later use.
if not os.path.exists(storage_dir):
    os.mkdir(storage_dir)
np.save('%s/ytree_redshifts' % (storage_dir), ytree_redshifts)
np.save('%s/ds_redshifts' % (storage_dir), ds_redshifts)
np.save('%s/cosmological_times' % (storage_dir), cosmological_times)


