#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import glob
import numpy as np
from netCDF4 import Dataset
import xarray as xr
import pickle
import os

# Notes:
# ======
# - Manual step: generate observational error data file "mean_error.nc"
# - Inversion steps at the end are not clear. Need new variable names. Possible to rewrite that
#   section in a more mathematically transparent way, or will that slow things down?
#
# - Removed variable observational error and replaced with fixed 15 ppb error. 
#   mean_error.nc (or mean_error_test.nc) is no longer needed.

# ==================================================================================================
#
#                                      Define functions
#
# ==================================================================================================

def save_obj(obj, name ):
    """ Save something with Pickle. """

    with open(name , 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


# --------------------------------------------------------------------------------------------------

def load_obj(name):
    """ Load something with Pickle. """

    with open( name, 'rb') as f:
        return pickle.load(f)


# --------------------------------------------------------------------------------------------------

def nearest_loc(loc_query, loc_grid, tolerance=1):
    """ Find the index of the nearest grid location to a query location, with some tolerance. """

    distances = np.abs(loc_grid - loc_query)
    ind = distances.argmin()
    if distances[ind] >= tolerance:
        return np.nan
    else:
        return ind


# --------------------------------------------------------------------------------------------------

#def get_prior_scaling(week, unit_sf_pth, cache):
    #sf = xr.load_dataset(unit_sf_pth)
    #for w in range(1,week):
    #    sf_w_pth = os.path.join(cache,f'inversion_result_week{w}.nc')
    #    sf_w = xr.load_dataset(sf_w_pth)
    #    sf['SF_Nonwetland'] = sf['SF_Nonwetland'] * sf_w['SF_Nonwetland']
    #prior_scaling = 1/sf['SF_Nonwetland'].values  # This seems wrong -- need to index by cluster id
    #return prior_scaling

def get_prior_scaling_constant(week, cache, repetitions, nclust=243):
    '''
    We need to multiply each row of the Jacobian by scalars = prior_new/prior_old (a vector with length nclust).
    For the constant prior case, prior_old is from the dynamic prior inversion and prior_new is the constant prior.
    So we have prior_old = a*prior_new, where a is a product of posterior scaling factors from different weeks in the dynamic prior inversion. 
    Thus scalars = 1/a
    
    Examples:
    - To get the scalars for week 5 so that we can solve the constant-prior inversion for week 5, input week here should be 5-1 = 4.
      The posteriors for weeks 1,2,3,4 will then be multiplied together, etc. 
    - The scalars for the week 1 (input week = 1-1 = 0) constant-prior inversion should be np.ones, because in week 1 the priors are identical in the 
      dynamic and constant cases.
    - The scalars for the week 2 (input week = 2-1 = 1) constant-prior inversion should come from the posterior scaling factors of the first week
      of the dynamic-prior inversion.
    '''
    if week == 0:
        prior_scaling_matrix = np.ones([repetitions,nclust])
    else:
        for w in range(1,week+1):
            if w == 1:
                sf = xr.open_dataset(f'{cache}/inversion_result_week{w}.nc')
            else:
                sf_w = xr.open_dataset(f'{cache}/inversion_result_week{w}.nc')
                sf['xhat'] = sf['xhat'] * sf_w['xhat']
        prior_scaling_vector = 1/sf['xhat'].values
        prior_scaling_matrix = np.tile(prior_scaling_vector, (repetitions,1))
   
    return prior_scaling_matrix
 

def get_prior_scaling_blended(week, cache_dynamic, cache_blended, repetitions, nclust=243, nbuff=8):
   '''
   Here we hold the buffer elements constant but allow the priors within the Permian to vary from week to week.
   Outside the Permian, the scalars are again = 1/a
   Inside the Permian, however, we have prior_old = a*edf_prior, and prior_new = b*edf_prior, 
   where edf_prior is the original prior and a and b are both products of posterior scaling factors:
   a is the product of scaling factors from the dynamic-prior inversion
   b is the product of scaling factors from the blended-prior inversion
   Hence prior_old = (a/b) * prior_new
   Hence within the Permian, the scalars are b/a
   '''
   if week == 0:
        prior_scaling_matrix = np.ones([repetitions,nclust])
   else:
       for w in range(1,week+1):
           if w == 1:
               sf_a = xr.open_dataset(f'{cache_dynamic}/inversion_result_week{w}.nc')
               sf_b = xr.open_dataset(f'{cache_blended}/inversion_result_week{w}.nc')
           else:
               sf_a_w = xr.open_dataset(f'{cache_dynamic}/inversion_result_week{w}.nc')
               sf_b_w = xr.open_dataset(f'{cache_blended}/inversion_result_week{w}.nc')
               sf_a['xhat'] = sf_a['xhat'] * sf_a_w['xhat']
               sf_b['xhat'] = sf_b['xhat'] * sf_b_w['xhat']
       # Reset scaling outside Permian to 1 for blended
       sf_b['xhat'][-nbuff:] = 1

       prior_scaling_vector = sf_b['xhat'].values / sf_a['xhat'].values
       
       prior_scaling_matrix = np.tile(prior_scaling_vector, (repetitions,1))

   return prior_scaling_matrix


def get_prior_scaling_uninformative(week, cache_dynamic, cache_uninformative, repetitions, nclust=243):
   '''
   This is the same as for the blended case, but without resetting the buffer element scaling factors to 1
   '''
   if week == 0:
        prior_scaling_matrix = np.ones([repetitions,nclust])
   else:
       for w in range(1,week+1):
           if w == 1:
               sf_a = xr.open_dataset(f'{cache_dynamic}/inversion_result_week{w}.nc')
               sf_b = xr.open_dataset(f'{cache_uninformative}/inversion_result_week{w}.nc')
           else:
               sf_a_w = xr.open_dataset(f'{cache_dynamic}/inversion_result_week{w}.nc')
               sf_b_w = xr.open_dataset(f'{cache_uninformative}/inversion_result_week{w}.nc')
               sf_a['xhat'] = sf_a['xhat'] * sf_a_w['xhat']
               sf_b['xhat'] = sf_b['xhat'] * sf_b_w['xhat']

       prior_scaling_vector = sf_b['xhat'].values / sf_a['xhat'].values

       prior_scaling_matrix = np.tile(prior_scaling_vector, (repetitions,1))

   return prior_scaling_matrix


def get_prior_scaling_nudged(week, cache_dynamic, cache_nudged, repetitions, nclust=243, nbuff=8):
    '''
    The same as for the blended case -- WRONG!! Need to emulate the nudging with original prior... 
    '''
    if week == 0:
        prior_scaling_matrix = np.ones([repetitions,nclust])
    else:
        for w in range(1,week+1):
            if w == 1:
                sf_a = xr.open_dataset(f'{cache_dynamic}/inversion_result_week{w}.nc')
                #sf_b = xr.open_dataset(f'{cache_nudged}/inversion_result_week{w}.nc')
            else:
                sf_a_w = xr.open_dataset(f'{cache_dynamic}/inversion_result_week{w}.nc')
                sf_a['xhat'] = sf_a['xhat'] * sf_a_w['xhat']
        post = xr.load_dataset(f'{cache_nudged}/archive_sf/sf_week{week}.nc') # WEEK + 1 ???? No!
        clusters = xr.load_dataset('/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/input_data_permian/Clusters_permian_kmeans.nc')
        sf_b = []
        for c in range(1,nclust+1):
            post_sf = np.nanmean(post['SF_Nonwetland'].where(clusters['Clusters'] == c).values)
            sf_b.append(post_sf)
        sf_b = np.asarray(sf_b)
        # Reset scaling outside Permian to 1 for nudged
        sf_b[-nbuff:] = 1

        prior_scaling_vector = sf_b / sf_a['xhat'].values

        prior_scaling_matrix = np.tile(prior_scaling_vector, (repetitions,1))

    return prior_scaling_matrix


# ==================================================================================================
#
#                                      Run the code
#
# ==================================================================================================

if __name__ == '__main__':
    import sys

    outputname = sys.argv[1]
    jacobian_dir = sys.argv[2]
    prior_opt = sys.argv[3] 
    week = int(sys.argv[4])    # not used for dynamic-prior inversion
    cache = sys.argv[5]        # not used for dynamic-prior inversion

    # Configuration
    n_clust = 235+8
    xlim = [-110,-96]#[-111,-95];
    ylim = [25,38]#[24,39]
    gamma = 0.3#0.4#0.25
    #jacobian_dir = "/n/holyscratch01/jacob_lab/dvaron/data_converted/"
    
    # Read observational error data
    #filename = "/n/holyscratch01/jacob_lab/dvaron/mean_error_test.nc"
    #data = xr.open_dataset(filename)
    #lon_GC = data['lon'].values
    #lat_GC = data['lat'].values
    #mean_error_std = data['error'].values
    #mean_error = mean_error_std**2                  # Error variance from standard deviation
    #mean_error = np.einsum('ij->ji', mean_error)
    #data.close()

    # Process the data_converted pkl files
    # ------------------------------------

    # Read output data from Step1 (virtual TROPOMI column XCH4, Jacobian matrix)
    os.chdir(jacobian_dir)
    files = glob.glob("*.pkl")
    files.sort()

    # Initialize ________
    all_part1 = np.zeros([n_clust,n_clust], dtype=float)
    all_part2 = np.zeros([n_clust], dtype=float)

    # For each .pkl file from Step1:
    for fi in files:
    
        # Load TROPOMI/GEOS-Chem and Jacobian matrix data from the .pkl file
        print(fi)
        met = load_obj(fi)
        # If there aren't any TROPOMI observations on this day, skip 
        # [****Shouldn't we have no files for those cases anyway?]
        if met['obs_GC'].shape[0] == 0:
            continue
        # Otherwise, grab the TROPOMI/GEOS-Chem data
        obs_GC = met['obs_GC']
        # Only consider data within latitude and longitude bounds
        ind = np.where((obs_GC[:,2]>=xlim[0]) & (obs_GC[:,2]<=xlim[1]) & (obs_GC[:,3]>=ylim[0]) & (obs_GC[:,3]<=ylim[1]))
        if (len(ind[0]) == 0):          # Skip if no data in bounds
            continue
        obs_GC = obs_GC[ind[0],:]       # TROPOMI and GEOS-Chem data within bounds
        KK = 1e9 * met['KK'][ind[0],:]  # Jacobian entries for observations within bounds [ppb]

        # If using dynamic prior
        if prior_opt == 'dynamic':
            # Do nothing
            pass

        # If using constant prior
        elif prior_opt == 'constant':
            repetitions = KK.shape[0]
            prior_scaling = get_prior_scaling_constant(week-1, cache, repetitions)
            KK = KK * prior_scaling

        # If using blended prior
        elif prior_opt == 'blended':
            blended_cache = '/'.join(outputname.split('/')[0:-1])
            repetitions = KK.shape[0]
            prior_scaling = get_prior_scaling_blended(week-1, cache, blended_cache, repetitions)
            KK = KK * prior_scaling

        # If using uninformative prior
        elif prior_opt == 'uninformative':
            uninformative_cache = '/'.join(outputname.split('/')[0:-1])
            repetitions = KK.shape[0]
            prior_scaling = get_prior_scaling_uninformative(week-1, cache, uninformative_cache, repetitions)
            KK = KK * prior_scaling

        # If using nudged prior
        elif prior_opt == 'nudged':
            nudged_cache = '/'.join(outputname.split('/')[0:-1])
            repetitions = KK.shape[0]
            prior_scaling = get_prior_scaling_nudged(week-1, cache, nudged_cache, repetitions)
            KK = KK * prior_scaling

        else:
            raise Exception("prior_opt must be dynamic, constant, blended, uninformative, or nudged.")
        NN = obs_GC.shape[0]            # Number of observations

        print('Sum of Jacobian entries:',np.sum(KK))

        # Now lower the sensitivity to BC by 50%
        # [****Is this something we want to do? Delete these lines?]
        #KK[:,1199:] = KK[:,1199:]/2
    
        # Initialize observation error matrix diagonal entries
        obs_error = np.zeros((NN,))
        # For each TROPOMI observation:
        for iNN in range(NN):
            ## Get the closest GC pixel by latitude and longitude
            #iGC = nearest_loc(obs_GC[iNN,2], lon_GC)
            #jGC = nearest_loc(obs_GC[iNN,3], lat_GC)
            ## Get the observational error for that location
            #obs_error[iNN] = mean_error[iGC,jGC]
            obs_error[iNN] = 15**2   # fixed error of 15 ppb
 
        # Measurement-model mismatch: TROPOMI columns minus GEOS-Chem virtual TROPOMI columns 
        deltaY = obs_GC[:,0] - obs_GC[:,1] # [ppb]
        # If there are any nan's in the data, abort 
        if (np.any(np.isnan(deltaY)) or np.any(np.isnan(KK)) or np.any(np.isnan(obs_error))):
            print('missing values', fi)
            break
    
        # [****What are these steps?]
        # [****We're clearly building some terms for the inversion here]
        KK_t = KK.transpose() 
        KK_t2 = np.zeros(KK_t.shape, dtype=float)
        for k in range(KK_t.shape[1]):
            KK_t2[:,k] = KK_t[:,k]/obs_error[k]        

        # [****What are part1 and part2?]
        # [****What is happening here?]
        part1 = KK_t2@KK
        part2 = KK_t2@deltaY
   
        # Add part1 & part2 to sums 
        all_part1 += part1
        all_part2 += part2
        
    # Prior covariance matrix
    # -----------------------

    # Relative errors, 1-sigma
    if prior_opt == 'uninformative':
        prior_err = 1e18
    else:
        prior_err = 0.5
    emis_error = np.zeros(n_clust);
    emis_error.fill(prior_err)      
 
    # Combined relative and absolute errors?
    combined_abs_rel_err = False
    if combined_abs_rel_err:
        print('Using combined absolute and relative errors')
        rel_err = 0.5
        abs_err = (4.9e-9)/10   # 10% of maximum prior emission
        clusters_path = '/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/input_data_permian/Clusters_permian_kmeans.nc'
        clusters = xr.load_dataset(clusters_path)
        prior_dir = '/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/production_permian/CH4_Jacobian/run_dirs/CH4_Jacobian_0000/OutputDir/'
        orig_prior_emis = xr.load_dataset(os.path.join(prior_dir,'HEMCO_diagnostics.201805010000.nc'))
        if prior_opt == 'blended':
            sf_path = '/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/sensitivity_permian/blended_prior_absrelerr/inversion/posterior_SF.nc'
            sf = xr.load_dataset(sf_path)
        else:
            # ERROR
            raise Exception("prior_opt must be blended.")
        current_prior = orig_prior_emis['EmisCH4_Total'].isel(time=0) * sf['SF_Nonwetland']
        abs_err_as_rel = np.abs(abs_err/current_prior)
        emis_error_array = rel_err + abs_err_as_rel
        emis_error = np.zeros(n_clust)
        for idx in range(1,n_clust+1):
            emis_error[idx-1] = np.nanmean(emis_error_array.where(clusters['Clusters'] == idx).values)
        emis_error[-8:] = 0.5

    # Build Sa matrix
    inv_Sa = np.diag(1/emis_error**2)   # Inverse of prior covariance matrix
    
    # Solve
    # -----
 
    # Solve for posterior scaling factors
    # [****Is that what "ratio" is?]
    ratio = np.linalg.inv(gamma*all_part1 + inv_Sa)@(gamma*all_part2)
    xhat = 1 + ratio # xhat = x_A + ratio = 1 + ratio

    # Posterior error
    #post_err = np.linalg.inv(all_part1 + inv_Sa)

    # Print some statistics
    print('Min:',xhat.min(),'Mean:',xhat.mean(),'Max',xhat.max())

    # Dictionary to store results
    # [****But we don't save it out? Delete these lines?]
    #met = {}
    #met['all_part1'] = all_part1
    #met['all_part2'] = all_part2
    #met['ratio'] = ratio

    # Save results
    #outputname = '/n/holyscratch01/jacob_lab/dvaron/inversion_result.nc'
    dataset = Dataset(outputname, 'w', format='NETCDF4_CLASSIC')
    nvar = dataset.createDimension('nvar', n_clust)
    nc_all_part1 = dataset.createVariable('all_part1', np.float32,('nvar','nvar'))
    nc_all_part2 = dataset.createVariable('all_part2', np.float32,('nvar'))
    nc_ratio = dataset.createVariable('ratio', np.float32,('nvar'))
    nc_xhat = dataset.createVariable('xhat', np.float32, ('nvar'))
    nc_all_part1[:,:] = all_part1    # [****What is this?]
    nc_all_part2[:] = all_part2      # [****What is this?]
    nc_ratio[:] = ratio              # [****What is this?]
    nc_xhat[:] = xhat
    dataset.close()
