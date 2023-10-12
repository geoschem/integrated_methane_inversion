# TODO: merge this script with invert.py to avoid redundancy
# This script performs the inversion but using lognormal 
# errors instead of normal errors. As an alternative to 
# the invert.py script. 

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os
import yaml
import pickle as pickle
import cartopy.crs as ccrs
import colorcet as cc
from scipy.sparse import spdiags
from datetime import datetime

# Implementing lognormal errors based on Zichong Chen's matlab script lrme1.m

# %cd /n/holylfs05/LABS/jacob_lab/shancock/imi_SA/inversion

config = yaml.load(open("/Users/lucasestrada/Downloads/Test_Permian_1week_14_0_2/config_Test_Permian_1week_14_0_2.yml"), Loader=yaml.FullLoader)
state_vector_filepath = "/Users/lucasestrada/Downloads/Test_Permian_1week_14_0_2/StateVector.nc"
inversion_data_pth = "/Users/lucasestrada/Downloads/Test_Permian_1week_14_0_2/inversion/"
state_vector = xr.load_dataset(state_vector_filepath)
state_vector_labels = state_vector['StateVector']
lon_bounds = [np.min(state_vector.lon.values), np.max(state_vector.lon.values)]
lat_bounds = [np.min(state_vector.lat.values), np.max(state_vector.lat.values)]
lats, lons = state_vector_labels.lat, state_vector_labels.lon

# GET TROPOMI DATA

ds = np.load(inversion_data_pth+"obs_ch4_tropomi.npz")
y_TROPOMI = np.asmatrix(ds["obs_ch4_tropomi"])
ds = np.load(inversion_data_pth+"gc_ch4_bkgd.npz")
ybkg_TROPOMI = np.asmatrix(ds["gc_ch4_bkgd"])
# ds = np.load(inversion_data_pth+"sve.npz")
# sve = ds["sve"]
y_ybkg_diff_TROPOMI = y_TROPOMI-ybkg_TROPOMI
# ds = np.load(inversion_data_pth+"albedo.npz", "wb")
# albedo = np.asmatrix(ds["albedo"])

# mask = (y_ybkg_diff_TROPOMI > 0) & (~np.isnan(sve)) & (albedo > 0.05)
# mask_simple = np.array(mask)[0]
# y_TROPOMI, ybkg_TROPOMI, y_ybkg_diff_TROPOMI = y_TROPOMI[mask], ybkg_TROPOMI[mask], y_ybkg_diff_TROPOMI[mask]

ds = np.load(inversion_data_pth+"K_blended.npz")
K600_TROPOMI = np.asmatrix(ds["K"])*1e9
ds = np.load(inversion_data_pth+"K16_blended.npz")
K16_TROPOMI = np.asmatrix(ds["K16"])*1e9

ds = np.load(inversion_data_pth+"so_super.npz")
so_TROPOMI = ds["so"]**2

# CONCATENATE GOSAT AND TROPOMI DATA
y = y_TROPOMI
ybkg = ybkg_TROPOMI
so = so_TROPOMI
K600 = K600_TROPOMI
K16 = K16_TROPOMI

y_ybkg_diff = y-ybkg
ybkg, y,  y_ybkg_diff = np.swapaxes(ybkg, 0,1), np.swapaxes(y, 0,1), np.swapaxes(y_ybkg_diff, 0,1)
K = np.concatenate((K600,K16),axis=1)

kappa=10
m,n = np.shape(K600)

xa = np.ones((n,1))*1.0
lnxa = np.asmatrix(np.log(xa))
xa_bcs = np.zeros((16,1))*1.0
xa = np.asmatrix(np.concatenate((xa,xa_bcs),axis=0))
lnxa = np.asmatrix(np.concatenate((lnxa,xa_bcs),axis=0))

soinv = 1/so
So = spdiags(so,0,m,m)
Soinv = spdiags(soinv,0,m,m)

lnsa_vals = [np.log(2), np.log(2.5), np.log(1.5)]
# lnsa_vals = [np.log(2)]
sa_bc_vals = [10,20,5]
# sa_bc_vals = [10]
# gamma_vals = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07]
# gamma_vals = [0.075,0.085,0.095]
gamma_vals = [0.025,0.02,0.3]
wetland_priors = ["unfccc","LPJ_CRU", "LPJ_ERA5", "LPJ_ERA5_MSWEP", "LPJ_MERRA2", "unfccc_wetcharts_2010_2019", ]
# wetland_priors = ["unfccc"]
for gamma in gamma_vals:
    invso_over_gamma = 1/(so/gamma)
    So_over_gamma_inv = spdiags(invso_over_gamma,0,m,m)
    for wetland_prior in wetland_priors:
        ds = np.load(inversion_data_pth+f"x_{wetland_prior}_blended.npz")
        xch40_TROPOMI =  np.asmatrix(ds["x_unfccc"][mask_simple])
        
        ds = np.load(inversion_data_pth+f"x_{wetland_prior}_gosat.npz")
        xch40_GOSAT =  np.asmatrix(ds["x_unfccc"])
        
        xch40 = np.append(xch40_GOSAT, xch40_TROPOMI, axis = 1)
        for lnsa_val in lnsa_vals:
            for sa_bc in sa_bc_vals:
                results_save_path = f"/n/holyscratch01/jacob_lab/shancock/results/results_gamma_{gamma}_lnsa_{lnsa_val}_bcsa_{sa_bc}_{wetland_prior}_with_gosat.npz"
                if os.path.exists(results_save_path):
                    continue
                print(gamma, lnsa_val, sa_bc, wetland_prior)
                lnsa = lnsa_val**2 * np.ones((n,1))
                sa_bcs = sa_bc**2*np.ones((16,1))
                lnsa_arr = np.concatenate((lnsa, sa_bcs), axis=0)
                # invlnsa = 1/lnsa
                # invlnsa_short = spdiags(np.squeeze(invlnsa[:600]),0,n,n)
                lnsa = np.zeros((n+16,n+16))
                np.fill_diagonal(lnsa, lnsa_arr)
                # lnsa = spdiags(np.squeeze(lnsa),0,n+16,n+16)
                # invlnsa = spdiags(np.squeeze(invlnsa),0,n+16,n+16)
                invlnsa = np.linalg.inv(lnsa)

                # INVERSION
                # Starting from prior
                xn=xa
                lnxn=lnxa

                lnk = np.concatenate((np.multiply(K600,np.transpose(xn[:600])), K16), axis = 1)

                print(len(y_ybkg_diff) - np.shape(y_ybkg_diff[y_ybkg_diff<0])[1])

                gamma_lnk_transpose_Soinv = gamma*np.transpose(lnk)@Soinv #Used in term1 and term2

                term1 = np.linalg.inv(gamma_lnk_transpose_Soinv@lnk+np.multiply((1+kappa),invlnsa))

                term2 = gamma_lnk_transpose_Soinv@(y_ybkg_diff-(K@xn))

                term3 = -1*invlnsa @ (lnxn - lnxa)

                ii = -1
                lnxn_update = lnxn+term1@(term2+term3)

                temp = max(abs(np.exp(lnxn_update[:600]) - np.exp(lnxn[:600]))/np.exp(lnxn[:600]))

                while temp >= 5e-3:
                    ii = ii+1
                    print("{}, The max relative diff is {}".format(ii, temp))
                    print("START: ", ii, datetime.now())
                    lnxn = lnxn_update
                    xn = np.concatenate((np.exp(lnxn[:600]), lnxn[600:]), axis = 0)
                    lnk = np.concatenate((np.multiply(K600,np.transpose(xn[:600])), K16), axis = 1)
                    gamma_lnk_transpose_Soinv = gamma*np.transpose(lnk)@Soinv #
                    term1 = np.linalg.inv(gamma_lnk_transpose_Soinv@lnk+(1+kappa)*invlnsa)
                    term2 = gamma_lnk_transpose_Soinv@(y_ybkg_diff-(K@xn))
                    term3 = -1*invlnsa @ (lnxn - lnxa)
                    lnxn_update = lnxn+term1@(term2+term3)
                    temp = max(abs(np.exp(lnxn_update[:600]) - np.exp(lnxn[:600]))/np.exp(lnxn[:600]))

                print("Done Iterating")
                lnxn=lnxn_update
                xn = np.concatenate((np.exp(lnxn[:600]), lnxn[600:]),axis=0)
                print("Post Proccessing")
                lnk = np.concatenate((np.multiply(K600,np.transpose(xn[:600])), K16), axis = 1)
                kso= np.transpose(lnk)@So_over_gamma_inv
                lns=np.linalg.inv(kso@lnk+invlnsa)
                G=lns@kso
                ak=G@lnk
                dofs=np.trace(ak);

                dlns=np.diag(lns[:600,:600]);
                xnmean= np.concatenate((np.multiply(xn[:600],np.expand_dims(np.exp(dlns*(0.5)), axis=1)),xn[600:]))

                Ja = np.transpose(lnxn[:600]-lnxa[:600])@invlnsa[:600,:600]@(lnxn[:600]-lnxa[:600])

                print("JA:", Ja)
                Ja_with_BCs = np.transpose(lnxn-lnxa)@invlnsa@(lnxn-lnxa)
                print("Ja with BCs", Ja_with_BCs)
                print(np.amin(xnmean[:600]),np.average(xnmean[:600]),np.amax(xnmean[:600]))

                # Plot up results
                xhat = xnmean
                xhat_arr = np.zeros((len(lats), len(lons)))
                for i in range(np.shape(xhat)[0]):
                    idx = np.where(state_vector_labels == float(i+1))
                    xhat_arr[idx] = xhat[i]
                xhat_arr = xr.DataArray(data=xhat_arr, 
                     dims=["lat", "lon"],
                     coords=dict(
                         lon=(["lon"], lons.values),
                         lat=(["lat"], lats.values),
                     ),)
                scale = xhat_arr
                fig = plt.figure(figsize=(8,8))
                plt.rcParams.update({'font.size': 16})
                ax = fig.subplots(1,1,subplot_kw={'projection': ccrs.PlateCarree()})

                plot_field(ax, scale, cmap='RdBu_r',
                           lon_bounds=lon_bounds, lat_bounds=lat_bounds,
                           vmin=-0.5, vmax=2.5, title='Scale factors', cbar_label='Scale factor',
                           only_ROI=False, state_vector_labels=state_vector_labels)

                ak_sensitivities = np.diagonal(ak)
                ak_arr = np.zeros((len(lats), len(lons)))
                for i in range(np.shape(ak_sensitivities)[0]):
                    idx = np.where(state_vector_labels == float(i+1))
                    ak_arr[idx] = ak_sensitivities[i]

                ak_arr = xr.DataArray(data=ak_arr, 
                                         dims=["lat", "lon"],
                                         coords=dict(
                                             lon=(["lon"], lons.values),
                                             lat=(["lat"], lats.values),
                                         ),)
                fig = plt.figure(figsize=(8,8))
                plt.rcParams.update({'font.size': 16})
                ax = fig.subplots(1,1,subplot_kw={'projection': ccrs.PlateCarree()})

                plot_field(ax, ak_arr, cmap=cc.cm.CET_L19,
                           lon_bounds=lon_bounds, lat_bounds=lat_bounds,
                           title='Averaging kernel sensitivities', cbar_label='Sensitivity', 
                           only_ROI=False, state_vector_labels=state_vector_labels)
                plt.show()

                np.savez(results_save_path, xn=xnmean,lnxn=lnxn,lns=lns,ak=ak,dofs=dofs,Ja=Ja)

                # Save out gridded posterior
                ds = xr.Dataset(
                    {
                        "ScaleFactor": (["lat", "lon"], scale.data)
                    },
                    coords={"lon": ("lon", lons.data), "lat": ("lat", lats.data)},
                )

                # Add attribute metadata
                ds.lat.attrs["units"] = "degrees_north"
                ds.lat.attrs["long_name"] = "Latitude"
                ds.lon.attrs["units"] = "degrees_east"
                ds.lon.attrs["long_name"] = "Longitude"
                ds.ScaleFactor.attrs["units"] = "1"

                # Create netcdf
                ds.to_netcdf(f"/n/holyscratch01/jacob_lab/shancock/results/gridded_posterior_gamma_{gamma}_lnsa_{lnsa_val}_bcsa_{sa_bc}_{wetland_prior}_with_gosat.nc")