# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 16:31:37 2024

@author: tilloal
"""
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import rasterio
import os, sys, inspect
import geopandas
import netCDF4 as nc
import scipy
import warnings
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")


main_path = 'D:/tilloal/Documents/06_Floodrivers/' ### CHANGE THIS PATH
valid_path = main_path + 'DataPaper/'
dis_path = main_path + 'HERA/'
os.chdir(valid_path)
EFAS_UpArea_dataset = rasterio.open(valid_path + 'GIS/upArea_European_01min.nc') ### EFAS upstream area file
EFAS_UpArea = EFAS_UpArea_dataset.read(1) / 1000000
Q_stations = geopandas.read_file('Europe_daily_combined_2022_v2_WGS84_NUTS.shp') #### SHP with gauge locations
Station_IDs = Q_stations.StationID
Countries = np.unique(Q_stations.Country.values)
years = list(range(1950, 2021))

#load the csv on matches to avoid spatial match on these stations
stat_efasmatch = pd.read_csv('Stations/efas_fl_match.csv').values
efm_ID=stat_efasmatch[:,0]


#load GRDC excel to match on upstream area when possible
stations_GRDC = pd.read_excel('GRDC_Stations.xlsx').values
grdc_ID=stations_GRDC[:,0].astype("int64")

# corr_all = np.full([Station_IDs.size,len(years)+2],np.inf)
# corr_all[:,0] = Station_IDs

# kge_all = np.full([Station_IDs.size,len(years)+2],np.inf)
# kge_all[:,0] = Station_IDs

# bias_all = np.full([Station_IDs.size,len(years)+2],np.inf)
# bias_all[:,0] = Station_IDs

# var_all = np.full([Station_IDs.size,len(years)+2],np.inf)
# var_all[:,0] = Station_IDs

# median_sim = np.full([Station_IDs.size,len(years)+2],np.inf)
# median_sim[:,0] = Station_IDs
# median_obs = np.full([Station_IDs.size,len(years)+2],np.inf)
# median_obs[:,0] = Station_IDs

# q1_sim = np.full([Station_IDs.size,len(years)+2],np.inf)
# q1_sim[:,0] = Station_IDs
# q1_obs = np.full([Station_IDs.size,len(years)+2],np.inf)
# q1_obs[:,0] = Station_IDs

# max_sim = np.full([Station_IDs.size,len(years)+2],np.inf)
# max_sim[:,0] = Station_IDs
# max_obs = np.full([Station_IDs.size,len(years)+2],np.inf)
# max_obs[:,0] = Station_IDs


# q2_sim = np.full([Station_IDs.size,len(years)+2],np.inf)
# q2_sim[:,0] = Station_IDs
# q2_obs = np.full([Station_IDs.size,len(years)+2],np.inf)
# q2_obs[:,0] = Station_IDs

#load the csv on matches to avoid spatial match on these stations
stat_corect = pd.read_csv('corrected_locations.csv').values
Station_IDs=stat_corect[:,0]
St_locs=stat_corect
for ky, y in enumerate(years):
    
    print(y)
    Q_data = pd.read_csv('Q_' + str(y) + '.csv', header=None).values  #### CSVs with observations
    Station_data_IDs = Q_data[0,]

    if y == 1950:
        offset = 3 # timesteps
        offsetX = 87 # remove the first three months of 1950
        offsetY = 90 
    else:
        offset = 1
        offsetX = 0
        offsetY = 0

    EFAS_Qann_dataset = rasterio.open('yearmean/Q_ymean'+str(y)+'.nc') ### average annual EFAS discharge files, see below
 	### USE THOSE BASH COMMANDS WITH CDO LOADED FIRST TO GENERATE ANNUAL AVERAGE DISCHARGE
 	#for i in {1950..2020}; do
 	#	cdo -z zip yearmean EFAS_outputs/Q_efas_final_${i}.nc EFAS_outputs/EFAS_${i}_avg.nc
 	#done

    EFAS_Qann = EFAS_Qann_dataset.read(1)

    EFAS_Qt_dataset = nc.Dataset(dis_path + 'dis.HERA_'+str(y)+'.nc') ### EFAS output files

    Dates = Q_data[1:, 0]

    times = int((EFAS_Qt_dataset.variables['time'].size)/4)
    Q_EFAS2 = np.full([times+1,Station_IDs.size],np.nan)
    Q_EFAS2[0, :] = Station_IDs
    T = range(0, times)
    for t in T:
        #I want to load all timesteps and then do a running average
        print(t)
        T_EFAS_t = EFAS_Qt_dataset['time'][t * 4+offset:(t + 1) * 4+offset]
        date_object = datetime(1979, 1, 1) + timedelta(days=T_EFAS_t[0])
        Q_EFAS_t = EFAS_Qt_dataset['dis'][t * 4+offset:(t + 1) * 4+offset, :, :]
        for k, s in enumerate(Station_IDs):
            X_EFAS = int(St_locs[k,1])-1
            Y_EFAS = int(St_locs[k,2])-1
            if X_EFAS > 0:
                suicide=Q_EFAS_t[:,Y_EFAS,X_EFAS]
                Q_EFAS2[t+1,k] = np.mean(Q_EFAS_t[:,Y_EFAS,X_EFAS])
    np.savetxt(valid_path + 'out/EFAS_corect_'+str(y)+'.csv', Q_EFAS2, fmt='%10.3f', delimiter=',')


#     Q_EFAS2=np.loadtxt(valid_path + 'out/EFAS_'+str(y)+'.csv', delimiter=',')
#     for k, s in enumerate(Station_IDs):
#         k = Station_IDs==s
#         ix = Station_data_IDs==s
#         Q_EFAS_s = np.squeeze(np.asarray(Q_EFAS2[1+offsetX:,k]))
#         if Q_EFAS_s.max() > 0:
#             Q_obs_s = np.squeeze(np.asarray(Q_data[1+offsetY:,ix][:,0]))
#             iy = np.isnan(Q_obs_s)
#             # Compute mean values excluding NaN
#             median_sim[k, ky+1] = np.nanmedian(Q_EFAS_s[~iy])
#             median_obs[k, ky+1] = np.nanmedian(Q_obs_s[~iy])
            
#             max_sim[k, ky+1] = np.nanmax(Q_EFAS_s[~iy])
#             max_obs[k, ky+1] = np.nanmax(Q_obs_s[~iy])
            
#             # if len(Q_obs_s[~iy])>0:
#             q1_sim[k, ky+1] = np.nanquantile(Q_EFAS_s[~iy],0.05)
#             q1_obs[k, ky+1] = np.nanquantile(Q_obs_s[~iy],0.05)
#             q2_sim[k, ky+1] = np.nanquantile(Q_EFAS_s[~iy],0.99)
#             q2_obs[k, ky+1] = np.nanquantile(Q_obs_s[~iy],0.99)
#             # else:
#             #     min_sim[k, ky+1] = np.nan
#             #     min_obs[k, ky+1] = np.nan
#             #     max_sim[k, ky+1] = np.nan
#             #     max_obs[k, ky+1] = np.nan
#             rp = scipy.stats.pearsonr(Q_EFAS_s[~iy], Q_obs_s[~iy])[0]
#             rs = scipy.stats.spearmanr(Q_EFAS_s, Q_obs_s, nan_policy='omit')[0]
#             b = np.mean(Q_EFAS_s[~iy])/np.mean(Q_obs_s[~iy])
#             g = (np.std(Q_EFAS_s[~iy])/np.mean(Q_EFAS_s[~iy]))/(np.std(Q_obs_s[~iy])/np.mean(Q_obs_s[~iy]))
#             bias_all[k, ky+1] = b
#             var_all[k, ky+1] =g
#             r2 = rs ** 2 if rs>0 else -(rs ** 2)
#             if np.isnan(rp):
#                 print("crap")
#                 corr_all[k, ky+1] = np.inf
#                 kge_all[k, ky+1] = 11
#             else:
#                 corr_all[k, ky+1] = rp
#                 kge_all[k, ky+1] = 1 - np.sqrt((rp-1)**2 + (b-1)**2 + (g-1)**2)
#             # Plot the array
#             # plt.plot(Q_EFAS_s, marker='o', linestyle='-')
#             # plt.plot(Q_obs_s, marker='o', linestyle='-')
#             # # Add labels and title
#             # plt.xlabel('Index')
#             # plt.ylabel('Values')
#             # plt.title('Array Plot')

# # # Show the plot
# # plt.show()

# corrs_countries = np.zeros([len(Countries)+1,len(years)+1])
# kge_countries = np.zeros([len(Countries)+1,len(years)+1])
# for k,c in enumerate(Countries):
#     Q_ix = Q_stations['Country']==c
#     c_c = corr_all[Q_ix.values, 1:len(years)+1]
#     k_c = kge_all[Q_ix.values, 1:len(years)+1]
#     ic = c_c != 0
#     ic2 = np.sum(ic,axis=1)>=0
#     corrs_countries[k, 0] = c_c[ic2,:].shape[0]
#     corrs_countries[k,1:] = np.mean(c_c[ic2,:],axis=0)
#     kge_countries[k, 0] = k_c[ic2,:].shape[0]
#     kge_countries[k,1:] = np.median(k_c[ic2,:],axis=0)

# icc = (kge_all!=0)&(kge_all<=1)
# kge_all[:,len(years)+1] = np.mean(kge_all,axis=1,where=icc)
# corr_all[:,len(years)+1] = np.mean(corr_all,axis=1,where=icc)

# np.savetxt(valid_path + 'out/EFAS_obs_KGE.csv', kge_all, fmt='%10.5f', delimiter=',')
# np.savetxt(valid_path + 'out/EFAS_obs_corr_pearson.csv', corr_all, fmt='%10.5f', delimiter=',')
# np.savetxt(valid_path + 'out/EFAS_obs_bias.csv', bias_all, fmt='%10.5f', delimiter=',')
# np.savetxt(valid_path + 'out/EFAS_obs_variability.csv', var_all, fmt='%10.5f', delimiter=',')

# np.savetxt(valid_path + 'out/EFAS_medianY.csv', median_sim, fmt='%10.5f', delimiter=',')
# np.savetxt(valid_path + 'out/obs_medianY.csv', median_obs, fmt='%10.5f', delimiter=',')

# np.savetxt(valid_path + 'out/EFAS_maxY.csv', max_sim, fmt='%10.5f', delimiter=',')
# np.savetxt(valid_path + 'out/obs_maxY.csv', max_obs, fmt='%10.5f', delimiter=',')

# np.savetxt(valid_path + 'out/EFAS_q05Y.csv', q1_sim, fmt='%10.5f', delimiter=',')
# np.savetxt(valid_path + 'out/obs_q05Y.csv', q1_obs, fmt='%10.5f', delimiter=',')
# np.savetxt(valid_path + 'out/EFAS_q99Y.csv', q2_sim, fmt='%10.5f', delimiter=',')
# np.savetxt(valid_path + 'out/obs_q99Y.csv', q2_obs, fmt='%10.5f', delimiter=',')

# np.savetxt(valid_path + 'out/EFAS_obs_correlations_countries_mean.csv', corrs_countries, fmt='%10.5f', delimiter=',')
# np.savetxt(valid_path + 'out/EFAS_obs_KGE_countries_mean.csv', kge_countries, fmt='%10.5f', delimiter=',')