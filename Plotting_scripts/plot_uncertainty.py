# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 16:25:42 2019

@author: jakob
"""
import numpy as np
import netCDF4 as nt
import matplotlib.pyplot as plt
from cartopy import config
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import xarray as xr
import time

start = time.time()

rootgrps = []
rootgrps.append(nt.Dataset("EN.4.2.1.analyses.g10.1950\EN.4.2.1.f.analysis.g10.195001.nc", "r+", format="NETCDF4"))
rootgrps.append(nt.Dataset("EN.4.2.1.analyses.g10.2018\EN.4.2.1.f.analysis.g10.201801.nc", "r+", format="NETCDF4"))

def moving_avg(times, HC,s):
    m = np.zeros((len(HC)-2*s))
    times = times[s:len(m)+s]
    for i in range(len(m)):
        n = i+s
        m[i] = np.mean(HC[n-s:n+s])
        
    return times, m
    
def calculate_weights(start_year, end_year, depth_from, depth_to):
    year = start_year
    month = 1 #(start month)
    years = [] #array to keep track of the years analysed (mainly for graph titles)
    times = [] #for the final plot of HC against time
    count = 0
    total_HC = np.zeros(((end_year-start_year+1)*12))
    DS = xr.open_mfdataset("EN.4.2.1.analyses.g10.1950\EN.4.2.1.f.analysis.g10.195001.nc",combine='by_coords',decode_times=False)
    
    for i in range((end_year-year+1)):
        years.append(year)
        for j in range(12):
            times.append(year+j/12)
            if month < 10: #months 1-9 are specified as "01"-"09" in the file names
                DS = xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+"0"+str(month)+".nc",combine='by_coords',decode_times=False)
                weights = DS.temperature_observation_weights
                HC_grid = weights.sel(depth=slice(depth_from,depth_to)).mean(dim="depth", skipna=True)
                total_HC[count] = HC_grid.isel(time=0).mean(skipna=True)

            else:
                DS = xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+str(month)+".nc",combine='by_coords',decode_times=False)
                weights = DS.temperature_observation_weights
                HC_grid = weights.sel(depth=slice(depth_from,depth_to)).mean(dim="depth", skipna=True)
                total_HC[count] = HC_grid.isel(time=0).mean(skipna=True)

            month += 1
            count += 1
        print("Year "+str(year)+" time elapsed = ", time.time() - start )
        year+=1
        month = 1 #reset to the first month of the year
    np.savetxt("HC_G_weights_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times,total_HC),delimiter=", ")
    times2, total_HC = moving_avg(times, total_HC, 12)
    np.savetxt("HC_G_MA_weights_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_HC),delimiter=", ")
    return times, total_HC, total_HC

calculate_weights(1950, 2018, 0, 700)
calculate_weights(1950, 2018, 700, 2000)


lon = rootgrps[0]["lon"][:]
lat = rootgrps[0]["lat"][:]
Tuncert = rootgrps[0]["temperature_uncertainty"][0,25,:,:]
weights = rootgrps[0]["temperature_observation_weights"][0,25,:,:]

#    clevs = np.linspace(-np.max(HC_slope), np.max(HC_slope), 50)   #contour levels

Tuncert_std = 8*np.std(Tuncert)
clevs = np.linspace(0, Tuncert_std, 20)   #contour levels
weights_std = 6*np.std(weights)
wclevs = np.linspace(0, weights_std, 20)   #contour levels

Tuncert, lon1 = add_cyclic_point(Tuncert, coord=lon)
weights, lon11 = add_cyclic_point(weights, coord=lon) 

Tuncert2 = rootgrps[1]["temperature_uncertainty"][0,25,:,:]
weights2 = rootgrps[1]["temperature_observation_weights"][0,25,:,:]

#    clevs = np.linspace(-np.max(HC_slope), np.max(HC_slope), 50)   #contour levels
Tuncert2_std = 8*np.std(Tuncert2)
clevs2 = np.linspace(0, Tuncert2_std, 20)   #contour levels
weights2_std = 6*np.std(weights2)
wclevs2 = np.linspace(0, weights2_std, 20)   #contour levels

Tuncert2, lon2 = add_cyclic_point(Tuncert2, coord=lon)
weights2, lon22 = add_cyclic_point(weights2, coord=lon) 


proj = ccrs.PlateCarree(180)
land_hires = cfeature.NaturalEarthFeature('physical', 'land', '50m', facecolor='grey')
fig = plt.figure(figsize=(10, 8))
ax1 = plt.subplot(2, 2, 1, projection=proj)
ax1.set_title("Temperature uncertainty in 1950 at depth = 1000m")
cf = plt.contourf(lon1, lat, Tuncert[:, :], levels=clevs, cmap=plt.cm.Reds, transform = ccrs.PlateCarree(0))
cb = plt.colorbar(cf, shrink=0.5, pad=0.02, orientation='vertical', fraction=0.1, format='%.1f')
cb.set_label('Std dev [°C]')
ax1.set_xticks([0, 60, 120, 180, 240, 300, 360], crs=ccrs.PlateCarree())
ax1.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)
ax1.coastlines()
ax1.add_feature(land_hires)
ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
              linewidth=1, color='gray', alpha=0.5, linestyle='--')


ax2 = plt.subplot(2, 2, 2, projection=proj)
ax2.set_title("Temperature uncertainty in 2018 at depth = 1000m")
cf2 = plt.contourf(lon2, lat, Tuncert2[:, :], levels=clevs2, cmap=plt.cm.Reds, transform = ccrs.PlateCarree(0))
cb2 = plt.colorbar(cf2, shrink=0.5, pad=0.02, orientation='vertical', fraction=0.1, format='%.1f')
cb2.set_label('Std dev [°C]')
ax2.coastlines()
ax2.add_feature(land_hires)
ax2.set_xticks([0, 60, 120, 180, 240, 300, 360], crs=ccrs.PlateCarree())
ax2.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax2.xaxis.set_major_formatter(lon_formatter)
ax2.yaxis.set_major_formatter(lat_formatter)
ax2.coastlines()
ax2.add_feature(land_hires)
ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
              linewidth=1, color='gray', alpha=0.5, linestyle='--')

ax3 = plt.subplot(2, 2, 3, projection=proj)
ax3.set_title("Observation weights 1950 at depth = 1000m")
cf3 = plt.contourf(lon11, lat, weights[:, :], cmap=plt.cm.Reds, transform = ccrs.PlateCarree(0))
cb3 = plt.colorbar(cf3, shrink=0.5, pad=0.02, orientation='vertical', fraction=0.1, format='%.1f')
cb3.set_label('Relative weighting')
ax3.coastlines()
ax3.add_feature(land_hires)
ax3.set_xticks([0, 60, 120, 180, 240, 300, 360], crs=ccrs.PlateCarree())
ax3.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax3.xaxis.set_major_formatter(lon_formatter)
ax3.yaxis.set_major_formatter(lat_formatter)
ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
              linewidth=1, color='gray', alpha=0.5, linestyle='--')

ax4 = plt.subplot(2, 2, 4, projection=proj)
ax4.set_title("Observation weights 2018 at depth = 1000m")
cf4 = plt.contourf(lon22, lat, weights2[:, :], cmap=plt.cm.Reds, transform = ccrs.PlateCarree(0))
cb4 = plt.colorbar(cf4, shrink=0.5, pad=0.02, orientation='vertical', fraction=0.1, format='%.1f')
cb4.set_label('Relative weighting')
ax4.coastlines()
ax4.add_feature(land_hires)
ax4.set_xticks([0, 60, 120, 180, 240, 300, 360], crs=ccrs.PlateCarree())
ax4.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax4.xaxis.set_major_formatter(lon_formatter)
ax4.yaxis.set_major_formatter(lat_formatter)
ax4.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
              linewidth=1, color='gray', alpha=0.5, linestyle='--')
plt.tight_layout()