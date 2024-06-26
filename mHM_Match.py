# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:05:43 2024

@author: tilloal
"""

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
dis_path = main_path + 'dis/calibrated/filtered/Histo/'
os.chdir(valid_path)
EFAS_UpArea_dataset = rasterio.open(valid_path + 'GIS/upArea_European_01min.nc') ### EFAS upstream area file
EFAS_UpArea = EFAS_UpArea_dataset.read(1) / 1000000
Q_stations = geopandas.read_file('Europe_daily_combined_2022_v2_WGS84_NUTS.shp') #### SHP with gauge locations
Station_IDs = Q_stations.StationID
Countries = np.unique(Q_stations.Country.values)
years = list(range(1950, 2021))

# Load the Excel file into a DataFrame
#%%
#load MmH point locations
mhm_path=valid_path+"Revisions/mHM_EU"
df = pd.read_excel(mhm_path+'/european_catchments_with_filename.xlsx')
mHM_loc=df.filename
# Print the first few rows of the DataFrame
print(df.head())
#%%
#match these locations with observed stations
#load the csv on matches to avoid spatial match on these stations
stat_efasmatch = pd.read_csv('Stations/efas_fl_match.csv').values
efm_ID=stat_efasmatch[:,0]


#%%
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

# for ky, y in enumerate(years):
    
#     print(y)
#     Q_data = pd.read_csv('Q_' + str(y) + '.csv', header=None).values  #### CSVs with observations
#     Station_data_IDs = Q_data[0,]

#     if y == 1950:
#         offset = 3 # timesteps
#         offsetX = 87 # remove the first three months of 1950
#         offsetY = 90 
#     else:
#         offset = 1
#         offsetX = 0
#         offsetY = 0

#     EFAS_Qann_dataset = rasterio.open('yearmean/Q_ymean'+str(y)+'.nc') ### average annual EFAS discharge files, see below
#  	### USE THOSE BASH COMMANDS WITH CDO LOADED FIRST TO GENERATE ANNUAL AVERAGE DISCHARGE
#  	#for i in {1950..2020}; do
#  	#	cdo -z zip yearmean EFAS_outputs/Q_efas_final_${i}.nc EFAS_outputs/EFAS_${i}_avg.nc
#  	#done

#     EFAS_Qann = EFAS_Qann_dataset.read(1)

#EFAS_Qt_dataset = nc.Dataset(dis_path + 'dis.filtered'+str(y)+'.nc') ### EFAS output files

#     Dates = Q_data[1:, 0]
St_locs = np.empty([mHM_loc.size, 3], dtype=object)
for k,s in enumerate(mHM_loc):
    #ix = Station_data_IDs==s
    #indices_boolean_indexing = np.nonzero(k)[0]
    print(k)
    St_locs[k,0]=s
    UpA_mHM = df.Area_given_km2[k]
    X = df.LON[k]
    Y = df.LAT[k]
    X_EFAS = int((X + 25.25) * 60)
    Y_EFAS = int((72.25 - Y) * 60)
    if X_EFAS<=4530 and Y_EFAS<=2970:
        ### search for nearest upstream area in surrounding grid cells
        UpArea = EFAS_UpArea[Y_EFAS-1:Y_EFAS+2,X_EFAS-1:X_EFAS+2]
        UpArea_l = UpArea.flatten()
        X_ind = np.array([0,1,2,0,1,2,0,1,2])
        Y_ind = np.array([0,0,0,1,1,1,2,2,2])
        UpA_diff = ((UpArea_l - UpA_mHM)**2)**(1/2)/UpA_mHM
        iyf = UpA_diff==min(UpA_diff)
        if UpA_diff[iyf][0] <= 5:
            St_locs[k,2] = X_EFAS - 1 + X_ind[iyf][0]
            St_locs[k,1] = Y_EFAS - 1 + Y_ind[iyf][0]
    else:
        St_locs[k,2] = np.nan
        St_locs[k,1] = np.nan

np.savetxt(mhm_path + '/mHM_locations_EFAS_grid.txt', St_locs, fmt='%s', delimiter=',',)

#%%
#Now I have to match the MHM points with validation points
#try with location and upstream area

#%%
main_path = 'D:/tilloal/Documents/06_Floodrivers/' ### CHANGE THIS PATH
valid_path = main_path + 'DataPaper/'
dis_path = main_path + 'HERA/'
os.chdir(valid_path)
EFAS_UpArea_dataset = rasterio.open(valid_path + 'GIS/upArea_European_01min.nc') ### EFAS upstream area file
EFAS_UpArea = EFAS_UpArea_dataset.read(1) / 1000000
Q_stations = geopandas.read_file('Europe_daily_combined_2022_v2_WGS84_NUTS.shp') #### SHP with gauge locations
Station_IDs = Q_stations.StationID
Countries = np.unique(Q_stations.Country.values)
years = list(range(1952, 2021))

#load the csv on matches to avoid spatial match on these stations
stat_efasmatch = pd.read_csv('Stations/efas_fl_match.csv').values
efm_ID=stat_efasmatch[:,0]

corr_mhm_all = np.full([Station_IDs.size,2],np.inf)
corr_mhm_all[:,0] = Station_IDs

kge_sqrt_all = np.full([Station_IDs.size,2],np.inf)
kge_sqrt_all[:,0] = Station_IDs

bias_sqrt_all = np.full([Station_IDs.size,2],np.inf)
bias_sqrt_all[:,0] = Station_IDs

var_sqrt_all = np.full([Station_IDs.size,2],np.inf)
var_sqrt_all[:,0] = Station_IDs

#I initiate with first year

EFAS_Qann_dataset = rasterio.open('yearmean/Q_ymean'+str(1951)+'.nc') ### average annual EFAS discharge files, see below
	### USE THOSE BASH COMMANDS WITH CDO LOADED FIRST TO GENERATE ANNUAL AVERAGE DISCHARGE
	#for i in {1950..2020}; do
	#	cdo -z zip yearmean EFAS_outputs/Q_efas_final_${i}.nc EFAS_outputs/EFAS_${i}_avg.nc
	#done

EFAS_Qann = EFAS_Qann_dataset.read(1)

EFAS_Qt_dataset = nc.Dataset(dis_path + 'dis.HERA_'+str(1951)+'.nc') ### EFAS output files


times = int((EFAS_Qt_dataset.variables['time'].size)/4)
Q_EFAS_ay=np.loadtxt(valid_path + 'out/EFAS_'+str(1951)+'.csv', delimiter=',')
Q_EFAS_ay=pd.DataFrame(Q_EFAS_ay)

Q_data_ay = pd.read_csv('Q_' + str(1951) + '.csv', header=None).values  #### CSVs with observation
Q_data_ay=pd.DataFrame(Q_data_ay)


for ky, y in enumerate(years):
    print(y)
    Q_data = pd.read_csv('Q_' + str(y) + '.csv', header=None).values  #### CSVs with observations
    Station_data_IDs = Q_data[0,]
    Q_datam=pd.DataFrame(Q_data)
    Q_datam=Q_datam.iloc[1:]
    print(len(Q_datam))
    Q_data_ay = pd.concat([Q_data_ay, Q_datam])


    Q_EFAS2=np.loadtxt(valid_path + 'out/EFAS_'+str(y)+'.csv', delimiter=',')
    Q_EFASm=pd.DataFrame(Q_EFAS2)
    Q_EFASm=Q_EFASm.iloc[1:]
    print(len(Q_EFASm))
    Q_EFAS_ay = pd.concat([Q_EFAS_ay, Q_EFASm])


#save Q_efas_ay

offsetX = 87 # remove the first three months of 1950
offsetY = 90

offsetX = 1 # remove the first three months of 1950
offsetY = 1

for k, s in enumerate(Station_IDs):
    print(k)
    s=Station_IDs[k]
    # k = Station_IDs==s
    ix = Station_data_IDs==s
    Q_EFAS_s = np.squeeze(np.asarray(Q_EFAS_ay)[1+offsetX:,k])
    print(np.nanmax(Q_EFAS_s))
    # plt.plot(Q_EFAS_s, marker='o', linestyle='--')
    # plt.xlabel('Index')
    # plt.ylabel('Q_efas_s')
    # plt.title('Q_efas_s Plot')
    # plt.show()
    if np.nanmax(Q_EFAS_s) > 0:
        Q_obs_s = np.squeeze(np.asarray(Q_data_ay)[1+offsetY:,ix][:,0])
        iy = np.isnan(Q_obs_s) & np.isnan(Q_EFAS_s)
        #iz= np.isnan(Q_EFAS_s[~iy])
        Q_EFAS_li=np.sqrt(Q_EFAS_s[~iy])
        
        Q_obs_li=np.sqrt(Q_obs_s[~iy])
        mask = np.isfinite(Q_EFAS_li) & np.isfinite(Q_obs_li)
        Q_EFAS_li=Q_EFAS_li[mask]
        Q_obs_li=Q_obs_li[mask]
        # Compute mean values excluding NaN
        # median_sim[k, 1] = np.nanmedian(Q_EFAS_s[~iy])
        # median_obs[k, 1] = np.nanmedian(Q_obs_s[~iy])
        # # if len(Q_obs_s[~iy])>0:
        # q1_sim[k, 1] = np.nanquantile(Q_EFAS_s[~iy],0.05)
        # q1_obs[k, 1] = np.nanquantile(Q_obs_s[~iy],0.05)
        # q2_sim[k, 1] = np.nanquantile(Q_EFAS_s[~iy],0.95)
        # q2_obs[k, 1] = np.nanquantile(Q_obs_s[~iy],0.95)
        
        
        # else:
        #     min_sim[k, ky+1] = np.nan
        #     min_obs[k, ky+1] = np.nan
        #     max_sim[k, ky+1] = np.nan
        #     max_obs[k, ky+1] = np.nan
        
        merd1=Q_EFAS_s[~iy]
        merd2=Q_obs_s[~iy]
        rp = scipy.stats.pearsonr((Q_EFAS_li), (Q_obs_li))[0]
        rs = scipy.stats.spearmanr(Q_EFAS_s, Q_obs_s, nan_policy='omit')[0]
        b = np.mean(Q_EFAS_li)/np.mean(Q_obs_li)
        g = (np.std(Q_EFAS_li)/np.mean(Q_EFAS_li))/(np.std(Q_obs_li)/np.mean(Q_obs_li))
        bias_sqrt_all[k, 1] = b
        var_sqrt_all[k, 1] =g
        r2 = rs ** 2 if rs>0 else -(rs ** 2)
        if np.isnan(rp):
            corr_sqrt_all[k, 1] = np.nan
            kge_sqrt_all[k, 1] = np.nan
        else:
            corr_sqrt_all[k, 1] = rp
            kge_sqrt_all[k, 1] = 1 - np.sqrt((rp-1)**2 + (b-1)**2 + (g-1)**2)
    else:
        print("olala")
        
        # Plot the array
        # plt.plot(Q_EFAS_s[~iy][1:30], marker='o', linestyle='--')
        # plt.plot(Q_obs_s[~iy][1:30], marker='o', linestyle='-')
        # # Add labels and title
        # plt.xlabel('Index')
        # plt.ylabel('Values')
        # plt.title('Array Plot')
        # plt.show()

# # # Show the plot


# corrs_countries = np.zeros([len(Countries)+1,len(years)+1])
# kge_countries = np.zeros([len(Countries)+1,len(years)+1])
# for k,c in enumerate(Countries):
# Q_ix = Q_stations['Country']==c
# c_c = corr_all[Q_ix.values, 1:len(years)+1]
# k_c = kge_all[Q_ix.values, 1:len(years)+1]
# ic = c_c != 0
# ic2 = np.sum(ic,axis=1)>=0
# corrs_countries[k, 0] = c_c[ic2,:].shape[0]
# corrs_countries[k,1:] = np.mean(c_c[ic2,:],axis=0)
# kge_countries[k, 0] = k_c[ic2,:].shape[0]
# kge_countries[k,1:] = np.median(k_c[ic2,:],axis=0)

# icc = (kge_all!=0)&(kge_all<=1)
# kge_all[:,len(years)+1] = np.mean(kge_all,axis=1,where=icc)
# corr_all[:,len(years)+1] = np.mean(corr_all,axis=1,where=icc)

np.savetxt(valid_path + 'out/EFAS_obs_sqrtKGEAY.csv', kge_sqrt_all, fmt='%10.5f', delimiter=',')
np.savetxt(valid_path + 'out/EFAS_obs_sqrtcorr_pearsonAY.csv', corr_sqrt_all, fmt='%10.5f', delimiter=',')
np.savetxt(valid_path + 'out/EFAS_obs_sqrtbiasAY.csv', bias_sqrt_all, fmt='%10.5f', delimiter=',')
np.savetxt(valid_path + 'out/EFAS_obs_sqrtvariabilityAY.csv', var_sqrt_all, fmt='%10.5f', delimiter=',')

# np.savetxt(valid_path + 'out/EFAS_obs_correlations_countries_mean.csv', corrs_countries, fmt='%10.5f', delimiter=',')
# np.savetxt(valid_path + 'out/EFAS_obs_KGE_countries_mean.csv', kge_countries, fmt='%10.5f', delimiter=',')


# times = int((EFAS_Qt_dataset.variables['time'].size)/4)
# Q_EFAS2 = np.full([times+1,Station_IDs.size],np.nan)
# Q_EFAS2[0, :] = Station_IDs
# T = range(0, times)
# for t in T:
#     #I want to load all timesteps and then do a running average
#     print(t)
#     T_EFAS_t = EFAS_Qt_dataset['time'][t * 4+offset:(t + 1) * 4+offset]
#     date_object = datetime(1979, 1, 1) + timedelta(days=T_EFAS_t[0])
#     Q_EFAS_t = EFAS_Qt_dataset['dis'][t * 4+offset:(t + 1) * 4+offset, :, :]
#     for k, s in enumerate(Station_IDs):
#         X_EFAS = int(St_locs[k,2])
#         Y_EFAS = int(St_locs[k,1])
#         if X_EFAS > 0:
#             suicide=Q_EFAS_t[:,Y_EFAS,X_EFAS]
#             Q_EFAS2[t+1,k] = np.mean(Q_EFAS_t[:,Y_EFAS,X_EFAS])
# np.savetxt(valid_path + 'out/EFAS_'+str(y)+'.csv', Q_EFAS2, fmt='%10.3f', delimiter=',')


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