import os
from os.path import join
import sys
import math

sys.path.append('.')
import inversion_lucas as inv

import xarray as xr
import numpy as np

import pandas as pd

######################
### FILE LOCATIONS ###
######################
main = '/Users/hannahnesser/Documents/Harvard/Research/Reduced_Rank_Jacobian'
code = join(main, 'python')
inputs = join(main, 'input')

#################
### LOAD DATA ###
#################
# Import clusters
clusters = xr.open_dataarray(join(inputs, 'clusters_1x125_plot.nc'))

# Load estimated and true Jacobian
k_est = xr.open_dataarray(join(inputs, 'k_est.nc'))
k_true = xr.open_dataarray(join(inputs, 'k_true.nc'))

# Load estimated dofs (you'll calculate this from K_est)
dofs = xr.open_dataarray(...)

# Define the native resolution state vector
nstate = len(dofs)
sv = np.arange(1, nstate + 1, 1)

#####################
### SET CONSTANTS ###
#####################
RF = 5
threshold = 0

########################################
### BUILD REDUCED DIMENSION JACOBIAN ###
########################################
# Generate 1 element state vector
sv = aggregate_cells(clusters, sv, dofs, 
                     n_cells=[nstate], n_cluster_size=[nstate])

# A record of the DOFS per cluster
dpc_record = [dofs.sum()/nstate]
dofs_record = [dofs.sum()]

# A record of the cluster size and the number of native-resolution grid
# cells that are used for that cluster size
cluster_size = [1]
cells = [25]

# A record of the state vector indices at which the cluster size increases
state_thresholds = [0]

# Initialize our while loop variables
dpc_check = True
size_count = 1
iteration_count = 1

# While the DOFS per cluster are less than the threshold set at the top
# of this script
while dpc_check:
    print('-'*100)
    print('Iteration : ', iteration_count)
    print('Max. Cluster Size: ', cluster_size[-1])

    # Update the state vector, construct the Jacobian, and solve the
    # inversion
    sv = inv.aggregate_cells(clusters, sv, dofs, n_cells=cells,
                             n_cluster_size=cluster_size)

    # Lucas: You will need to recalculate the DOFS for the new state
    # vector. You should be able to do this using the W matrix 
    # approach described in most of our group's inversion papers.
    dofs = update_dofs(dofs) # (This function doesn't actually exist)

    # Add the number of state vector elements and the DOFS per cluster
    # to the records
    nstate_record.append(np.sum(cells))
    dpc_record.append(dofs.sum()/sv.max())
    dofs_record.append(dofs.sum())

    # Check whether the change in DOFS per cluster from the previous
    # estimate is greater than the threshold
    dpc_check = ((dpc_record[-1] - dpc_record[-2]) > threshold)
                 # or (iteration_count <= 1))
    sv_check = (np.sum(cells) < est0.nstate)

    # If all of the state vector elements are allocated to the state
    # vector, stop the loop
    if not sv_check:
        print('ALL NATIVE-RESOLUTION STATE VECTOR ELEMENTS ASSIGNED.')
        print(f'ITERATIONS : {iteration_count}')
        state_thresholds.append(np.sum(cells))
        dpc_check = False
    # If it is less than the threshold, increase the cluster size
    elif not dpc_check:
        print('Updating cluster size.')
        # Increase the cluster size and add more native-resolution
        # state vector elements to be allocated to the state vector
        # Add to state_thresholds
        state_thresholds.append(np.sum(cells))
        cluster_size.append(2**size_count)
        cells.append(cluster_size[-1]*cell_base)

        # Update dpc_check
        dpc_check = True

        # Up the counts
        size_count += 1
    else:
        # Add more native-resolution state vector elements to be
        # allocated
        cells[-1] += cluster_size[-1]*cell_base

    # Up the iteration count
    iteration_count += 1

print('Complete.')
