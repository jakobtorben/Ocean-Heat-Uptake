# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:47:28 2019

@author: benro

NB: lat = lat index - 83
    lon = lon index + 1
"""

import numpy as np
import os
os.environ["PROJ_LIB"] = "C:/Users/benro/.conda/pkgs/proj4-5.2.0-ha925a31_1/Library/share"; #fixr
import time
import xarray as xr
import matplotlib.pyplot as plt

start = time.time()
Cp = 3850 #Heat capacity of seawater
rho = 1025 #density of seawater


def calculate_HC_global(start_year, end_year, depth_from, depth_to):
    year = start_year
    month = 1 #(start month)
    years = [] #array to keep track of the years analysed (mainly for graph titles)
    times = [] #for the final plot of HC against time
    count = 0
    total_HC = np.zeros(((end_year-start_year+1)*12))
    for i in range((end_year-year+1)):
        years.append(year)
        for j in range(12):
            times.append(year+j/12)
            if month < 10: #months 1-9 are specified as "01"-"09" in the file names
                #rootgrps.append(nt.Dataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+"0"+str(month)+".nc", "r+", format="NETCDF4"))
                DS = xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+"0"+str(month)+".nc",combine='by_coords',decode_times=False)
                DS.values
                DS.time.values
                temps = DS.temperature.fillna(0)
                depths = DS.sel(depth=slice(depth_from,depth_to))["depth"].values
                #temps.sel(lat=0, lon=slice(320,359), depth=slice(0,1000)).plot(hue = 'lon')
                temps_interp = temps.sel(depth=slice(depth_from,depth_to)).interp(depth = np.arange(depths[0],depths[-1],1), method = 'cubic')
                #HC_grid_i = temps_interp.sel(depth=slice(0,1000)).integrate("depth")*Cp*rho
                HC_grid_i = temps_interp.sel(depth=slice(depth_from,depth_to)).sum("depth")*Cp*rho
                #print(HC_grid_i)
                temps = temps.assign_coords(HC_grid_i = HC_grid_i)
                area = 6400000**2*(np.cos(temps.lat * np.pi/180)) * (temps.lon/temps.lon) * (np.pi/180)**2
                temps = temps.assign_coords(area = area)
                HC_grid_weighted = (temps.HC_grid_i * temps.area)
                total_HC[count] = HC_grid_weighted.isel(time=0).values.sum()
                print('Total heat content: ' + str(total_HC[count]))
            else:
                #rootgrps.append(nt.Dataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+str(month)+".nc", "r+", format="NETCDF4"))
                DS = xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+str(month)+".nc",combine='by_coords',decode_times=False)
                DS.values
                DS.values
                DS.time.values
                temps = DS.temperature.fillna(0)
                depths = DS.sel(depth=slice(depth_from,depth_to))["depth"].values
                #temps.sel(lat=0, lon=slice(320,359), depth=slice(0,1000)).plot(hue = 'lon')
                temps_interp = temps.sel(depth=slice(depth_from,depth_to)).interp(depth = np.arange(depths[0],depths[-1],1), method = 'cubic')
                #HC_grid_i = temps_interp.sel(depth=slice(0,1000)).integrate("depth")*Cp*rho
                HC_grid_i = temps_interp.sel(depth=slice(depth_from,depth_to)).sum("depth")*Cp*rho
                #print(HC_grid_i)
                temps = temps.assign_coords(HC_grid_i = HC_grid_i)
                area = 6400000**2*(np.cos(temps.lat * np.pi/180)) * (temps.lon/temps.lon) * (np.pi/180)**2
                temps = temps.assign_coords(area = area)
                HC_grid_weighted = (temps.HC_grid_i * temps.area)
                total_HC[count] = HC_grid_weighted.isel(time=0).values.sum()
                print('Total heat content: ' + str(total_HC[count]))
            month += 1
            count += 1
        year+=1
        
        month = 1 #reset to the first month of the year
        
    return times, total_HC


def monthly_avgs(HC):
    """For a given dataset of HCs, outputs an array of length 12 giving the mean HC value 
    for each month."""
    monthlies = np.zeros((int(len(HC)/12),12))
    counter_m = 0   #keeps track of years
    counter_n = 0   #keeps track of months
    for i in range(len(HC)):
        if counter_n<12:
            monthlies[counter_m,counter_n] = HC[i]
            counter_n += 1
        else:
            counter_m += 1
            monthlies[counter_m,0] = HC[i]
            counter_n = 1
    monthly_avgs = np.zeros((12))
    months = np.zeros((12))
    for i in range(12):
        monthly_avgs[i] = np.mean(monthlies[:,i])
        months[i] = i+1
        
    return months, monthly_avgs


def run_global(start_year, end_year, depth_from, depth_to):
    """Basically a main; runs all functions defined above."""
    #years, times, rootgrps = retrieve(start_year,end_year)
    #four_temps, depths, lats, lons = unpack_rootgrps(rootgrps)
    times, total_HC = calculate_HC_global(start_year, end_year, depth_from, depth_to)
    
    plt.plot(times, total_HC)
    
    #months, month_avgs = monthly_avgs(HC)
    
    return HC #, months, month_avgs

#years, times, HC, pos, months, month_avgs = run(1950, 2018, 25, 31)
HC = run_global(1950, 1952, 0, 500)

print(time.time()-start)