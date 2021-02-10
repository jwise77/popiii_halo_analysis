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

# Load datasets and arbor
datadir = '/home/jwise/runs/SG64-2020/SG64-L2'
ts = yt.load('%s/DD????/output_????' % (datadir))
a = ytree.load('%s/rockstar_halos/trees/tree_0_0_0.dat' % (datadir))

# Determine number of datasets.
N = len(ts)

# Loop through halos in the arbor to find the progenitor line with the longest history
# This will be the root halo whose progenitor line will be used as reference for redshifts.
maxN = -1
ref_halo = None
for root_halo in a:
    if root_halo['prog', 'id'].size > maxN:
        maxN = root_halo['prog', 'id'].size
        ref_halo = root_halo
print("Storing %d of %d redshifts that have halo data" % (maxN, N))
N = maxN

# There are slight discrepancies between redshifts stored in the datasets and
# redshifts in the saved ytree arbor so save those as separate arrays.
ytree_redshifts = np.zeros(N)
ds_redshifts = np.zeros(N)
cosmological_times = np.zeros(N)

# Loop through each dataset and halo in the progenitor line to save redshifts
# and cosmological time.
# Reverse order of progenitor line to go in increasing time.
for i, ds in enumerate(ts[-N:]):
    ytree_redshifts[i] = ref_halo['prog', 'redshift'][-i]
    ds_redshifts[i] = ds.current_redshift
    cosmological_times[i] = ds.current_time.to('yr')

# Save all arrays to directory for later use.
if not os.path.exists('stored_arrays'):
    os.mkdir('stored_arrays')
np.save('stored_arrays/ytree_redshifts', ytree_redshifts)
np.save('stored_arrays/ds_redshifts', ds_redshifts)
np.save('stored_arrays/cosmological_times', cosmological_times)


