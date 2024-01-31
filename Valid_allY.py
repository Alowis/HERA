# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:30:00 2023

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
years = list(range(1951, 2021))

#load the csv on matches to avoid spatial match on these stations
stat_efasmatch = pd.read_csv('Stations/efas_fl_match.csv').values
efm_ID=stat_efasmatch[:,0]


#Initialisation of output matrices
corr_all = np.full([Station_IDs.size,2],np.inf)
corr_all[:,0] = Station_IDs

kge_all = np.full([Station_IDs.size,2],np.inf)
kge_all[:,0] = Station_IDs

bias_all = np.full([Station_IDs.size,2],np.inf)
bias_all[:,0] = Station_IDs

var_all = np.full([Station_IDs.size,2],np.inf)
var_all[:,0] = Station_IDs

mean_sim = np.full([Station_IDs.size,2],np.inf)
mean_sim[:,0] = Station_IDs
mean_obs = np.full([Station_IDs.size,2],np.inf)
mean_obs[:,0] = Station_IDs

q1_sim = np.full([Station_IDs.size,2],np.inf)
q1_sim[:,0] = Station_IDs
q1_obs = np.full([Station_IDs.size,2],np.inf)
q1_obs[:,0] = Station_IDs

q2_sim = np.full([Station_IDs.size,2],np.inf)
q2_sim[:,0] = Station_IDs
q2_obs = np.full([Station_IDs.size,2],np.inf)
q2_obs[:,0] = Station_IDs


percentiles = range(1,100)
qperf_sim = np.full([Station_IDs.size,100],np.inf)
qperf_sim[:,0] = Station_IDs
qperf_obs = np.full([Station_IDs.size,100],np.inf)
qperf_obs[:,0] = Station_IDs

#I initiate with first year
Q_EFAS_ay=np.loadtxt(valid_path + 'out/EFAS_'+str(1950)+'.csv', delimiter=',')
Q_EFAS_ay=pd.DataFrame(Q_EFAS_ay)
Q_data_ay = pd.read_csv('Q_' + str(1950) + '.csv', header=None).values  #### CSVs with observation
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


offsetX = 87 # remove the first three months of 1950
offsetY = 90

for k, s in enumerate(Station_IDs):
    print(k)
    ix = Station_data_IDs==s
    Q_EFAS_s = np.squeeze(np.asarray(Q_EFAS_ay)[1+offsetX:,k])
    print(np.nanmax(Q_EFAS_s))
    if np.nanmax(Q_EFAS_s) > 0:
        Q_obs_s = np.squeeze(np.asarray(Q_data_ay)[1+offsetY:,ix][:,0])
        iy = np.isnan(Q_obs_s)
        iz= np.isnan(Q_EFAS_s[~iy])
        Q_EFAS_li=Q_EFAS_s[~iy]
        Q_obs_li=Q_obs_s[~iy]
        
        # Compute mean values excluding NaN
        mean_sim[k, 1] = np.nanmean(Q_EFAS_li[~iz])
        mean_obs[k, 1] = np.nanmean(Q_obs_li[~iz])
        qperf_sim[k, np.arange(1, 100)] = np.nanpercentile(Q_EFAS_li[~iz],percentiles)
        qperf_obs[k, np.arange(1, 100)] = np.nanpercentile(Q_obs_li[~iz],percentiles)

        out1=Q_EFAS_s[~iy]
        out2=Q_obs_s[~iy]
        rp = scipy.stats.pearsonr(Q_EFAS_li[~iz], Q_obs_li[~iz])[0]
        rs = scipy.stats.spearmanr(Q_EFAS_s, Q_obs_s, nan_policy='omit')[0]
        b = np.mean(Q_EFAS_li[~iz])/np.mean(Q_obs_li[~iz])
        g = (np.std(Q_EFAS_li[~iz])/np.mean(Q_EFAS_li[~iz]))/(np.std(Q_obs_li[~iz])/np.mean(Q_obs_li[~iz]))
        bias_all[k, 1] = b
        var_all[k, 1] =g
        r2 = rs ** 2 if rs>0 else -(rs ** 2)
        if np.isnan(rp):
            corr_all[k, 1] = np.nan
            kge_all[k, 1] = np.nan
        else:
            corr_all[k, 1] = rp
            kge_all[k, 1] = 1 - np.sqrt((rp-1)**2 + (b-1)**2 + (g-1)**2)
    else:
        print("olala")
        mean_sim[k, 1] = np.nan
        mean_obs[k, 1] = np.nan
        
#Saving the outputs
np.savetxt(valid_path + 'out/EFAS_obs_KGEAY.csv', kge_all, fmt='%10.5f', delimiter=',')
np.savetxt(valid_path + 'out/EFAS_obs_corr_pearsonAY.csv', corr_all, fmt='%10.5f', delimiter=',')
np.savetxt(valid_path + 'out/EFAS_obs_biasAY.csv', bias_all, fmt='%10.5f', delimiter=',')
np.savetxt(valid_path + 'out/EFAS_obs_variabilityAY.csv', var_all, fmt='%10.5f', delimiter=',')
np.savetxt(valid_path + 'out/EFAS_meanAY.csv', mean_sim, fmt='%10.5f', delimiter=',')
np.savetxt(valid_path + 'out/obs_meanAY.csv', mean_obs, fmt='%10.5f', delimiter=',')
np.savetxt(valid_path + 'out/EFAS_percentilesAY.csv', qperf_sim, fmt='%10.5f', delimiter=',')
np.savetxt(valid_path + 'out/obs_percentilesAY.csv', qperf_obs, fmt='%10.5f', delimiter=',')

