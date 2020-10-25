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


def calculate_HC_global_ranges(start_year, end_year, start_depth, end_depth, savetxt):
    year = start_year
    month = 1 #(start month)
    years = [] #array to keep track of the years analysed (mainly for graph titles)
    times = [] #for the final plot of HC against time
    count = 0
    re = 6378.137*1e3 #equatorial radius earth (m)
    e = 0.08181919  #eccentricity
    degree = np.pi/180
    total_HC = np.zeros(((end_year-start_year+1)*12))
    
    if end_depth<=200:
        increment = 10
    elif end_depth<=300:
        increment = 30
    elif end_depth<=1000:
        increment = 50
        
    end_depth_sel = end_depth
    
    for i in range((end_year-year+1)):
        years.append(year)
        for j in range(12):
            times.append(year+j/12)
            if month < 10: #months 1-9 are specified as "01"-"09" in the file names
                DS = xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+"0"+str(month)+".nc",combine='by_coords',decode_times=False)
            else:
                DS = xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+str(month)+".nc",combine='by_coords',decode_times=False)
            DS.values
            DS.time.values
            temps = DS.temperature.fillna(0)
            depths = DS.sel(depth=slice(start_depth,end_depth_sel))["depth"].values
            while depths[-1]<end_depth:
                end_depth_sel += 1
                depths = DS.sel(depth=slice(start_depth,end_depth_sel))["depth"].values
            #temps.sel(lat=0, lon=slice(320,359), depth=slice(0,1000)).plot(hue = 'lon')
            temps_interp = temps.sel(depth=slice(start_depth,end_depth_sel)).interp(depth = np.arange(depths[0],depths[-1],1), method = 'cubic')
            
            depth_from = start_depth
            depth_to = start_depth+increment
            while depth_to <= end_depth:
                rtop = re-depth_from
                HC_grid_i = temps_interp.sel(depth=slice(depth_from,depth_to)).sum("depth")*Cp*rho
                
                temps = temps.assign_coords(HC_grid_i = HC_grid_i)
                area = rtop*rtop*(np.cos(temps.lat*degree))*(1-e*e)*(temps.lon/temps.lon)*degree*degree/(1-e*e*np.sin(temps.lat*degree)*np.sin(temps.lat*degree))**2
                temps = temps.assign_coords(area = area)
                HC_grid_weighted = (temps.HC_grid_i * temps.area)
                if savetxt == True:
                    np.savetxt("HC_grid_"+"year_"+str(year)+"_month_"+str(month)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",HC_grid_weighted.values[0,:,:],delimiter=", ")
            
                depth_from += increment
                depth_to += increment
                
                #total_HC[count] = HC_grid_weighted.isel(time=0).values.sum()
                #print('Total heat content: ' + str(total_HC[count]))
            print("    Month "+str(month)+" done")
            month += 1
            count += 1
        print("Year: "+str(year-start_year+1)+". Time elapsed: "+str(time.time()-start))
        year+=1
        month = 1 #reset to the first month of the year
        
    return times, total_HC


def moving_avg(times, HC,s):
    m = np.zeros((len(HC)-2*s))
    times = times[s:len(m)+s]
    for i in range(len(m)):
        n = i+s
        m[i] = np.mean(HC[n-s:n+s])
        
    return times, m


def run_global(start_year, end_year, start_depth, end_depth, savetxt = False):
    """Basically a main; runs all functions defined above."""
    #years, times, rootgrps = retrieve(start_year,end_year)
    #four_temps, depths, lats, lons = unpack_rootgrps(rootgrps)
    times, total_HC = calculate_HC_global_ranges(start_year, end_year, start_depth, end_depth, savetxt = savetxt)
    
    #times, total_HC = moving_avg(times, total_HC, 12)
    
    #total_HC = total_HC - np.mean(total_HC[:629])
    #plt.plot(times, total_HC)
    #plt.show()    
    return "done" #, months, month_avgs

#years, times, HC, pos, months, month_avgs = run(1950, 2018, 25, 31)

#HC = run_global(1950, 2018, 0, 200, savetxt=True)
HC = run_global(1950, 2018, 0, 1000, savetxt=True)

print(time.time()-start)