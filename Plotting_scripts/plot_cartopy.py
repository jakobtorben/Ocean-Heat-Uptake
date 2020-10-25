# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 10:54:04 2019

@author: jakob
"""
import numpy as np
from netCDF4 import Dataset
import os
import matplotlib.pyplot as plt
from cartopy import config
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point


d_1960_01 = Dataset('EN.4.2.1.f.analysis.g10.196001.nc')
depth = d_1960_01.variables['depth'][:]
lat = d_1960_01.variables['lat'][:]
lon = d_1960_01.variables['lon'][:]
temperature = d_1960_01.variables['temperature'][:]
temperature_uncertainty = d_1960_01.variables['temperature_uncertainty'][:]

clevs = np.linspace(np.min(temperature[0,0,:,:]),
                    np.max(temperature[0, 0,:,:]), 20)   #contour levels

temperature, lon = add_cyclic_point(temperature, coord=lon)

fig = plt.figure(figsize=(10,7))
proj = ccrs.PlateCarree()
ax = plt.axes(projection=proj)
cf = plt.contourf(lon, lat, temperature[0,0,:,:], levels=clevs, cmap=plt.cm.coolwarm, transform=proj)
cb = plt.colorbar(cf, extend='both', shrink=0.675, pad=0.02, orientation='vertical', fraction=0.1)
cb.ax.set_ylabel('Potential Temperature [K]')    
ax.coastlines()
ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.tick_params(direction='out')
plt.show()

