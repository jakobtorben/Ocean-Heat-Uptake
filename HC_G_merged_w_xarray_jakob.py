# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:47:28 2019

@author: benro

NB: lat = lat index - 83
    lon = lon index + 1
"""

import numpy as np
import time
import xarray as xr
import matplotlib.pyplot as plt

start = time.time()
Cp = 3850 #Heat capacity of seawater
rho = 1025 #density of seawater


def calculate_HC_global(start_year, end_year, depth_from, depth_to, uncertainty=False):
    year = start_year
    month = 1 #(start month)
    years = [] #array to keep track of the years analysed (mainly for graph titles)
    times = [] #for the final plot of HC against time
    count = 0
    total_HC = np.zeros(((end_year-start_year+1)*12))
    DS = xr.open_mfdataset("EN.4.2.1.analyses.g10.1950/EN.4.2.1.f.analysis.g10.195001.nc",combine='by_coords',decode_times=False)
    DS = DS.sel(depth=slice(depth_from,depth_to))
    temps = DS.temperature
    depths = DS["depth"].values
        
    
    re = 6378.137*1e3 #equatorial radius earth (m)
    rtop = re-depths[0]
    e = 0.08181919  #eccentricity
    dln = dlt = 1*np.pi/180
    lt = temps.lat*np.pi/180
    area = rtop*rtop*np.cos(lt)*(1-e*e)*dlt*dln/(1-e*e*np.sin(lt)*np.sin(lt))**2*(temps.lon/temps.lon)
    
    for i in range((end_year-year+1)):
        years.append(year)
        for j in range(12):
            times.append(year+j/12)
            if month < 10: #months 1-9 are specified as "01"-"09" in the file names
                DS = xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"/EN.4.2.1.f.analysis.g10."+str(year)+"0"+str(month)+".nc",combine='by_coords',decode_times=False)
                DS = DS.sel(depth=slice(depth_from, depth_to))
                temps = DS.temperature.fillna(0)
                temps_interp = temps.interp(depth = np.arange(depths[0],depths[-1],1), method = 'cubic')
                HC_grid_i = temps_interp.sum("depth")*Cp*rho
                temps = temps.assign_coords(HC_grid_i = HC_grid_i)
                temps = temps.assign_coords(area = area)
                HC_grid_weighted = (temps.HC_grid_i * temps.area)
                total_HC[count] = HC_grid_weighted.isel(time=0).values.sum()

            else:
                DS = xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"/EN.4.2.1.f.analysis.g10."+str(year)+str(month)+".nc",combine='by_coords',decode_times=False)
                DS = DS.sel(depth=slice(depth_from, depth_to))
                temps = DS.temperature.fillna(0)
                temps_interp = temps.interp(depth = np.arange(depths[0],depths[-1],1), method = 'cubic')
                HC_grid_i = temps_interp.sum("depth")*Cp*rho
                temps = temps.assign_coords(HC_grid_i = HC_grid_i)
                temps = temps.assign_coords(area = area)
                HC_grid_weighted = (temps.HC_grid_i * temps.area)
                total_HC[count] = HC_grid_weighted.isel(time=0).values.sum()
                
            month += 1
            count += 1
        print("Year "+str(year)+" time elapsed = ", time.time() - start )
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

def std_baseline(ttl_gbl):
    std_dev = np.std(ttl_gbl)
    mean = np.mean(ttl_gbl)
    base_yrs_end = np.where(ttl_gbl > (mean+std_dev))[0][0]
    ttl_gbl -= np.mean(ttl_gbl[:base_yrs_end])
    print(base_yrs_end)
    
    return ttl_gbl


def run_global(start_year, end_year, depth_from, depth_to):
    """Basically a main; runs all functions defined above."""
    times, total_HC = calculate_HC_global(start_year, end_year, depth_from, depth_to)

    #plt.plot(times, total_HC)
    
#    np.savetxt("HC_G_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times,total_HC),delimiter=", ")
#    np.savetxt("HC_G_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times,std_baseline(total_HC)),delimiter=", ")
#    times2, total_HC = moving_avg(times, total_HC, 12)
#    plt.plot(times2, total_HC)
#    np.savetxt("HC_G_MA_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2, total_HC),delimiter=", ")
#    total_HC = std_baseline(total_HC)
#    np.savetxt("HC_G_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_HC),delimiter=", ")

HC = run_global(1950, 1951, 0, 500)

print(time.time()-start)