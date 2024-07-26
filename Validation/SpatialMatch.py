# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 14:42:06 2023

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

#load the csv on matches to avoid spatial match on these stations
stat_efasmatch = pd.read_csv('Stations/efas_fl_match.csv').values
efm_ID=stat_efasmatch[:,0]
corr_all = np.zeros([Station_IDs.size,len(years)+2])
corr_all[:,0] = Station_IDs

kge_all = np.zeros([Station_IDs.size,len(years)+2])
kge_all[:,0] = Station_IDs

mean_sim = np.zeros([Station_IDs.size,len(years)+2])
mean_sim[:,0] = Station_IDs
mean_obs = np.zeros([Station_IDs.size,len(years)+2])
mean_obs[:,0] = Station_IDs

min_sim = np.zeros([Station_IDs.size,len(years)+2])
min_sim[:,0] = Station_IDs
min_obs = np.zeros([Station_IDs.size,len(years)+2])
min_obs[:,0] = Station_IDs

max_sim = np.zeros([Station_IDs.size,len(years)+2])
max_sim[:,0] = Station_IDs
max_obs = np.zeros([Station_IDs.size,len(years)+2])
max_obs[:,0] = Station_IDs

all_out=np.zeros([Station_IDs.size,6])

#ky and y
ky=1
y=1951
Q_data = pd.read_csv('Q_' + str(y) + '.csv', header=None).values  #### CSVs with observations
Station_data_IDs = Q_data[0,]

if y == 1950:
    offset = 3 # timesteps
    offset1 = 3 # days
    offsetX = 91 # remove the first three months of 1950
else:
    offset = 1
    offset1 = 0
    offsetX = 1

EFAS_Qann_dataset = rasterio.open('yearmean/Q_ymean'+str(y)+'.nc') ### average annual EFAS discharge files, see below
EFAS_Qann = EFAS_Qann_dataset.read(1)
EFAS_Qt_dataset = nc.Dataset(dis_path + 'dis.filtered'+str(y)+'.nc') ### EFAS output files
Dates = Q_data[1:, 0]

St_locs = np.zeros([Station_IDs.size, 3])
for k,s in enumerate(Station_IDs):
    ix = Station_data_IDs==s
    print(s)
    efm=efm_ID==s
    efas_meta=stat_efasmatch[efm,]
    St_locs[k,0]=s
    if len(efas_meta)>0:
        St_locs[k,2] = efas_meta[0][2]
        St_locs[k,1] = efas_meta[0][1]
    else:
        Q_s = Q_data[1:,ix]
        if len(Q_s)>0:
            iyd = np.isnan(Q_s)
            d = len(Q_s)-sum(iyd)
            if d>183:
                Q_s_a = np.mean(Q_s[~iyd])
                if Q_s_a > 0:
                    X = Q_stations.centroid[k].x
                    Y = Q_stations.centroid[k].y
                    X_EFAS = int((X + 25.25) * 60)
                    Y_EFAS = int((72.25 - Y) * 60)
                    ### search for nearest average annual discharge in surrounding grid cells
                    Q_EFAS = EFAS_Qann[Y_EFAS-1:Y_EFAS+2,X_EFAS-1:X_EFAS+2]
                    Q_EFAS_l = Q_EFAS.flatten()
                    UpArea = EFAS_UpArea[Y_EFAS-1:Y_EFAS+2,X_EFAS-1:X_EFAS+2]
                    UpArea_l = UpArea.flatten()
                    X_ind = np.array([0,1,2,0,1,2,0,1,2])
                    Y_ind = np.array([0,0,0,1,1,1,2,2,2])
                    Q_diff = ((Q_EFAS_l - Q_s_a)**2)**(1/2)/Q_s_a
                    iyf = Q_diff==min(Q_diff)
                    if UpArea_l[iyf][0] >= 100:
                        St_locs[k,2] = X_EFAS - 1 + X_ind[iyf][0]
                        St_locs[k,1] = Y_EFAS - 1 + Y_ind[iyf][0]
    np.savetxt(valid_path + 'Stations_locations_EFAS_grid_'+str(y)+'.txt', St_locs, fmt='%d', delimiter=',',)