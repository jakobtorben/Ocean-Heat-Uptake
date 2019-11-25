# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:47:28 2019

@author: benro

NB: lat = lat index - 83
    lon = lon index + 1
"""

import matplotlib.pyplot as plt
import numpy as np
import os
os.environ["PROJ_LIB"] = "C:/Users/benro/.conda/pkgs/proj4-5.2.0-ha925a31_1/Library/share"; #fixr
import time
from scipy.interpolate import interp1d
import netCDF4 as nt
import pandas as pd

import xarray as xr
import dask

start = time.time()
Cp = 3850 #Heat capacity of seawater
rho = 1025 #density of seawater


def retrieve(start_year,end_year):
    rootgrps = []
    year = start_year
    month = 1 #(start month)
    years = [] #array to keep track of the years analysed (mainly for graph titles)
    times = [] #for the final plot of HC against time
    for i in range((end_year-year+1)):
        years.append(year)
        for j in range(12):
            times.append(year+j/12)
            if month < 10: #months 1-9 are specified as "01"-"09" in the file names
                #rootgrps.append(nt.Dataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+"0"+str(month)+".nc", "r+", format="NETCDF4"))
                rootgrps.append(xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+"0"+str(month)+".nc",combine='by_coords',decode_times=False))
                rootgrps[len(rootgrps)-1].values
                rootgrps[len(rootgrps)-1].time.values
            else:
                #rootgrps.append(nt.Dataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+str(month)+".nc", "r+", format="NETCDF4"))
                rootgrps.append(xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+str(month)+".nc",combine='by_coords',decode_times=False))
                rootgrps[len(rootgrps)-1].values
                rootgrps[len(rootgrps)-1].time.values
            month += 1
        year+=1
        month = 1 #reset to the first month of the year
    return years, times, rootgrps

def unpack_rootgrps(rootgrps):
    n_months = int(len(rootgrps))
    n_depths = int(len(rootgrps[0]["temperature"][0,:,0,0]))
    n_lats = len(rootgrps[0]["temperature"][0,0,:,0])
    n_lons = len(rootgrps[0]["temperature"][0,0,0,:])
    
    four_temps = np.zeros((n_months, n_depths, n_lats, n_lons))
    depths = np.zeros((n_depths))
    lats = np.zeros((n_lats))
    lons = np.zeros((n_lons))
    for month in range(n_months):
        four_temps[month, :, :, :] = rootgrps[month]["temperature"][0,:,:,:]
    depths = rootgrps[0]["depth"]
    lats = rootgrps[0]["lat"]
    lons = rootgrps[0]["lon"]
        
    return four_temps, depths, lats, lons
            

def depth_ind(depths, depth_from, depth_to):
    """Input depth range in (m), outputs index as array with [0]=depth_from, [1]=depth_to"""
    count = 0
    for i in range(len(depths)):
        if depths[i] >= depth_from and count == 0:
            depth_from_i = i
            count += 1
        elif depths[i] > depth_to and depths[i-1] <= depth_to:
            depth_to_i = i - 1
    return depth_from_i, depth_to_i


def calculate_HC_global(start_year, end_year):
    year = start_year
    month = 1 #(start month)
    years = [] #array to keep track of the years analysed (mainly for graph titles)
    times = [] #for the final plot of HC against time
    for i in range((end_year-year+1)):
        years.append(year)
        for j in range(1):
            times.append(year+j/12)
            if month < 10: #months 1-9 are specified as "01"-"09" in the file names
                #rootgrps.append(nt.Dataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+"0"+str(month)+".nc", "r+", format="NETCDF4"))
                DS = xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+"0"+str(month)+".nc",combine='by_coords',decode_times=False)
                DS.values
                DS.time.values
                temps = DS.temperature.fillna(0)
                #temps_interp = temps.sel(depth=slice(0,1000)).interp(np.linspace(0,1000,1), method="cubic")
                #depth = temps.sel(lat=0,lon=1)
                #temps.sel(lat=0, lon=slice(320,359), depth=slice(0,1000)).plot(hue = 'lon')
                temps_interp = temps.sel(depth=slice(0,1000)).interp(method= 'cubic')
                temps_interp.isel(lat=83, lon=slice(320,321)).plot(hue = 'lon')
                HC_grid_i = temps_interp.sel(depth=slice(0,1000)).integrate("depth")*Cp*rho
                #print(HC_grid_i)
                temps = temps.assign_coords(HC_grid_i = HC_grid_i)
                area = 6400000**2 * (np.cos(temps.lat * np.pi/180)) * (temps.lon/temps.lon) * np.pi/180**2
                temps = temps.assign_coords(area = area)
                HC_grid_weighted = (temps.HC_grid_i * temps.area)
                print('Total heat content: ' + str(HC_grid_weighted.isel(time=0).values.sum()))
            else:
                #rootgrps.append(nt.Dataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+str(month)+".nc", "r+", format="NETCDF4"))
                DS = xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+str(month)+".nc",combine='by_coords',decode_times=False)
                DS.values
                DS.time.values
                temps = DS.temperature.fillna(0)
            month += 1
        year+=1
        month = 1 #reset to the first month of the year 
#    HC_grid_i = np.ma.zeros([len(four_temps[:,0,0,0]),len(lats),len(lons)])
#    HC_grid = np.ma.masked_all_like(HC_grid_i)
#    year = start_year
#    month = 1
#    for i in range(1):  #loops through every month
#        three_temps = four_temps[i,:,:,:]
#        names = ["x","y","z"]
#        index = pd.MultiIndex.from_product([range(s) for s in three_temps.shape], names = names)
#        df = pd.DataFrame({"Temps": three_temps.flatten()}, index = index)["Temps"]
#        
#        df = df.reorder_levels(["x","z","y"]).sort_index()
#        df = df.unstack(level = "x").swaplevel().sort_index()
#        df.columns = [str(i) for i in range(42)]
#        df.index.names = ["LAT", "LON"]
#        df = df.mask(df<=0)
#        
#        time_a = time.time()
#        
#        for row in df.itertuples():
#            row*=2
#        for lat in range(len(three_temps[:,0,0])):
#            for lon in range(len(three_temps[0,:,0])):
#                three_temps[lat,lon,:]*=2
            
#        print(time.time()-time_a)
        #print(df)
        
#    A = np.random.randint(0,1000,(3,5,4))
#    
#    names = ['x', 'y', 'z']
#    index = pd.MultiIndex.from_product([range(s)for s in A.shape], names=names)
#    df = pd.DataFrame({'A': A.flatten()}, index=index)['A']
#    
#    df = df.unstack(level='x').swaplevel().sort_index()
#    df.columns = ['A', 'B', 'C']
#    df.index.names = ['DATE', 'i']
#    
#    print(df)
        
#        for lat in range(len(lats)):
#            for lon in range(len(lons)):
#                temp_mask = four_temps[i,depth_from:(depth_to+1),lat,lon].mask
#                if np.ma.all(temp_mask) == False: #checks if all temp data at location are empty
#                    if np.ma.any(temp_mask) == True:   # finds data with empty temps to separate values depth<depthrange
#                        pos = np.where(temp_mask == True)[0][0]
#                        depths = depths_i[depth_from:(pos)] #finds the limited depth
#                    else:
#                        depths = depths_i[depth_from:(depth_to)] #depth is given by depthrange  
#                    temps = four_temps[i,depth_from:(depth_from+len(depths)),lat,lon]
#                    if len(temps) > 3:  #can only interpolate with enough data points
#                        cs = interp1d(depths,temps,kind='cubic')
#                        xs = np.arange(depths[0],depths[len(depths)-1],1) #depths spaced by 1m
#                        theta = cs(xs)
#                        ones = np.ones((len(xs)))  
#                        #HC_grid[lt, ln] = np.dot(theta,ones) #sums over all interpolated temperatures more efficiently
#                        HC_grid[i,lat,lon] = np.dot(theta,ones)
#        HC_i = HC_grid[i,:,:]*rho*Cp
#        np.savetxt("HC_grid_"+"year_"+str(year)+"_month_"+str(month)+"_"+str(depth_from_i)+"m_"+str(depth_to_i)+"m.txt",HC_i,delimiter=", ")
#        if month < 12:
#            month += 1
#        else:
#            month = 1
#            year += 1
#    HC_grid *= rho*Cp
#    return HC_grid


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


def run_global(start_year, end_year, depth_from, depth_to, animate=True):
    """Basically a main; runs all functions defined above."""
    #years, times, rootgrps = retrieve(1950,2018)
    #rootgrps = [nt.Dataset("EN.4.2.1.analyses.g10.1950\EN.4.2.1.f.analysis.g10.195001.nc", "r+", format="NETCDF4")]
    #years, times, rootgrps = retrieve(start_year,end_year)
    #four_temps, depths, lats, lons = unpack_rootgrps(rootgrps)
    HC = calculate_HC_global(start_year, end_year)
    
    #months, month_avgs = monthly_avgs(HC)
    
    return HC #, months, month_avgs

#years, times, HC, pos, months, month_avgs = run(1950, 2018, 25, 31)
HC = run_global(1950, 1950, 0, 500, False)

#fig, ax = plt.subplots(figsize=(6.5, 4))
#ax.plot(depths,temps0,"x")
#ax.plot(xs0,cs0(xs0))
#plt.title(str(years[0]))
#plt.grid()
#plt.show()
print(time.time()-start)