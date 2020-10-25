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
import os
import matplotlib.pyplot as plt
from cartopy.mpl.ticker import LatitudeFormatter
import matplotlib.ticker as ticker
from cartopy import config
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import cartopy.feature as cfeature


start = time.time()
Cp = 3850 #Heat capacity of seawater
rho = 1025 #density of seawater


def calculate_HC_global(start_year, end_year, depth_from, depth_to):
    year = start_year
    month = 1 #(start month)
    years = [] #array to keep track of the years analysed (mainly for graph titles)
    times = [] #for the final plot of HC against time
    count = 0
    lat = np.zeros(((end_year-start_year+1)*12, 173))
    DS = xr.open_mfdataset("EN.4.2.1.analyses.g10.1950\EN.4.2.1.f.analysis.g10.195001.nc",combine='by_coords',decode_times=False)
    depths = DS.sel(depth=slice(depth_from,depth_to))["depth"].values
    print(depths)
    temps = DS.temperature.fillna(0)
    
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
                DS = xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+"0"+str(month)+".nc",combine='by_coords',decode_times=False)
                temps = DS.temperature.fillna(0)
                temps_interp = temps.sel(depth=slice(depth_from,depth_to)).interp(depth = np.arange(depths[0],depths[-1],1), method = 'cubic')
                HC_grid_i = temps_interp.sel(depth=slice(depth_from,depth_to)).sum("depth")*Cp*rho
                temps = temps.assign_coords(HC_grid_i = HC_grid_i)
                temps = temps.assign_coords(area = area)
                HC_grid_weighted = (temps.HC_grid_i * temps.area)
                lat[count, :] = HC_grid_weighted.isel(time=0).sum(dim=('lon'))
            
            else:
                DS = xr.open_mfdataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+str(month)+".nc",combine='by_coords',decode_times=False)
                temps = DS.temperature.fillna(0)
                temps_interp = temps.sel(depth=slice(depth_from,depth_to)).interp(depth = np.arange(depths[0],depths[-1],1), method = 'cubic')
                HC_grid_i = temps_interp.sel(depth=slice(depth_from,depth_to)).sum("depth")*Cp*rho
                temps = temps.assign_coords(HC_grid_i = HC_grid_i)
                temps = temps.assign_coords(area = area)
                HC_grid_weighted = (temps.HC_grid_i * temps.area)
                lat[count, :] = HC_grid_weighted.isel(time=0).sum(dim=('lon'))

                
            month += 1
            count += 1
        
        print("Year "+str(year)+" time elapsed = ", time.time() - start )
        year+=1
        month = 1 #reset to the first month of the year
        
    latsarr = xr.DataArray(lat, coords=[times, range(-83,90)], dims=['time', 'lat'])
    latsarr = latsarr.mean(dim=('time'))
    return latsarr

def calculate_HC_global_ranges(start_year, end_year, start_depth, end_depth):
    year = start_year
    month = 1 #(start month)
    years = [] #array to keep track of the years analysed (mainly for graph titles)
    times = [] #for the final plot of HC against time
    count = 0
    re = 6378.137*1e3 #equatorial radius earth (m)
    e = 0.08181919  #eccentricity
    degree = np.pi/180
    increment = 50
    lat = np.zeros(((end_year-start_year+1)*12, 173, end_depth//increment))
    
#    if end_depth<=200:
#        increment = 50
#    elif end_depth<=300:
#        increment = 50
#    elif end_depth<=1000:
#        increment = 50
        
    end_depth_sel = end_depth
    
    for i in range((end_year-year+1)):
        years.append(year)
        for j in range(12):
            times.append(year+j/12)
            depth_count = 0
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
                lat[count, :, depth_count] = HC_grid_weighted.isel(time=0).sum(dim=('lon'))
                #if savetxt == True:
                 #   np.savetxt("HC_grid_"+"year_"+str(year)+"_month_"+str(month)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",HC_grid_weighted.values[0,:,:],delimiter=", ")
                depth_count += 1
                depth_from += increment
                depth_to += increment
                
            print("    Month "+str(month)+" done")
            month += 1
            count += 1
        print("Year: "+str(year-start_year+1)+". Time elapsed: "+str(time.time()-start))
        year+=1
        month = 1 #reset to the first month of the year
        
    latsarr = xr.DataArray(lat, coords=[times, range(-83,90), range(increment//2, int(end_depth+increment//2), increment)], dims=['time','lat', 'depth'])
    latsarr = latsarr.mean(dim=('time'))
    return latsarr


def run_global():
    """Basically a main; runs all functions defined above."""
    lat1 = calculate_HC_global_ranges(2004, 2008, 0, 2000)
    lat2 = calculate_HC_global_ranges(2014, 2018, 0, 2000)
    data = lat2 - lat1
    data.to_netcdf('zonal_2014-2018vs1950-1954_0to2000m.nc')

def plot(file): 

    ds = xr.open_dataset(str(file))
    data = ds.to_array() 
    
    lats = range(-83, 90)
    depths= data['depth'].values
    clevs = np.linspace(-5e20, 5e20, 41)
    
    @ticker.FuncFormatter
    def major_formatter(x, pos):
        label = str(-x) if x < 0 else str(x)
        return label


    fig = plt.figure(figsize=(11, 7))   
    ax = plt.subplot()
    depths, lats = np.meshgrid(depths, lats)
    cf = plt.contourf(lats, np.negative(depths), data[0,:,:], levels=clevs, cmap=plt.cm.bwr)
    cb = plt.colorbar(cf, shrink=0.9, pad=0.1, orientation='horizontal', fraction=0.1)
    cb.set_label('OHC anom. (J)')
    ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_yticks([0, -250, -500, -750, -1000, -1250, -1500, -1750, -2000])
    ax.yaxis.set_major_formatter(major_formatter)
    ax.set_ylabel("Depth (m)")
    ax.set_xlabel('Latitude')
    labelsx = [item.get_text() for item in ax.get_xticklabels()]
    labelsx[0] = u'90°S'
    labelsx[1] = u'60°S'
    labelsx[2] = u'30°S'
    labelsx[3] = u'0°'
    labelsx[4] = u'30°N'
    labelsx[5] = u'60°N'
    labelsx[6] = u'90°N'
    ax.set_xticklabels(labelsx)
    #fig.savefig("Figures\Linear_OHC_change_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.png")
    
    

#run_global()
plot('zonal_2014-2018vs1950-1954_0to2000m.nc')
plot('zonal_2014-2018vs2004-2008_0to2000m.nc')

print(time.time()-start)