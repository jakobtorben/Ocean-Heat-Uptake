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
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import interp1d
import netCDF4 as nt
#from time import sleep
import sys

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
                rootgrps.append(nt.Dataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+"0"+str(month)+".nc", "r+", format="NETCDF4"))
            else:
                rootgrps.append(nt.Dataset("EN.4.2.1.analyses.g10."+str(year)+"\EN.4.2.1.f.analysis.g10."+str(year)+str(month)+".nc", "r+", format="NETCDF4"))
            month += 1
        year+=1
        month = 1 #reset to the first month of the year
    return years, times, rootgrps



def depth_ind(rootgrps, depth_from, depth_to):
    """Input depth range in (m), outputs index as array with [0]=dpeth_from, [1]=depth_to"""
    depths = rootgrps[0]["depth"][:]
    count = 0
    for i in range(len(depths)):
        if depths[i] >= depth_from and count == 0:
            depth_from_i = i
            count += 1
        elif depths[i] > depth_to and depths[i-1] <= depth_to:
            depth_to_i = i - 1
    return depth_from_i, depth_to_i            


def calculate_HC(rootgrps,depth_from,depth_to,lat,lon,calculate_errors = False):
    depth_from, depth_to = depth_ind(rootgrps, depth_from, depth_to)
    depths = rootgrps[0]["depth"][depth_from:depth_to]
    xs = np.arange(depths[0],depths[len(depths)-1],1) #depths spaced by 1m
    ones = np.zeros((len(xs)))
    for i in range(len(ones)):
        ones[i] = 1
    HC = np.zeros((len(rootgrps)))
    
    if calculate_errors == False or depth_to>25:
        for i in range(len(rootgrps)):
            temps = rootgrps[i]["temperature"][0,depth_from:depth_to,lat+83,lon-1]
            cs = interp1d(depths,temps,kind='cubic')
            theta = cs(xs)
            HC[i] = np.dot(theta,ones) #sums over all interpolated temperatures more efficiently
        HC *= rho*Cp
        return HC
    
    else:
        HC_Err = np.zeros((len(rootgrps)))
        for i in range(len(rootgrps)):
            temps = rootgrps[i]["temperature"][0,depth_from:depth_to,lat+83,lon-1]
            temps_Err = rootgrps[i]["temperature_uncertainty"][0,depth_from:depth_to,lat+83,lon-1]
            cs = interp1d(depths,temps,kind='cubic')
            cs_err = interp1d(depths,temps_Err,kind='cubic')
            theta = cs(xs)
            #theta_err = cs_err(xs)
            HC[i] = np.dot(theta,ones)
            #HC_Err[i] = np.dot(theta_e)
        HC *= rho*Cp
        HC_Err *= rho*Cp
        return HC


def calculate_HC_global(rootgrps,depth_from,depth_to):
    depth_from, depth_to = depth_ind(rootgrps, depth_from, depth_to)
    lon = rootgrps[0]["lon"][:]
    lat = rootgrps[0]["lat"][:]

    ln2, lt2 = np.meshgrid(lon,lat)
    HC_grid_i = (np.ma.zeros(np.shape(lt2)))
    HC_grid = np.ma.masked_all_like(HC_grid_i)

    for i in range(len(rootgrps)):  #loops through every month
        for lt in range(len(lat)):#range(int(lat[0]), int(lat[len(lat)-1])+1):  #latitude
            for ln in range(len(lon)):#range(int(lon[0]), int(lon[len(lon)-1])+1):  #longitutde
                if np.ma.all(rootgrps[i]["temperature"][0,depth_from:(depth_to+1),lt,ln].mask) == False: #checks if all temp data at location are empty  
                    if np.ma.any(rootgrps[i]["temperature"][0,depth_from:(depth_to+1),lt,ln].mask) == True:   # finds data with empty temps to separate values depth<depthrange
                        count = 0
                        for j in rootgrps[i]["temperature"][0,depth_from:(depth_to+1),lt,ln].mask:
                            if j == True:
                                count += 1
                        depths = rootgrps[i]["depth"][depth_from:(depth_to-count+1)] #finds the limited depth
                    else:
                        depths = rootgrps[i]["depth"][depth_from:(depth_to+1)] #depth is given by depthrange  
                    temps = (rootgrps[i]["temperature"][0,depth_from:(depth_to+len(depths)),lt,ln])
                    if len(temps) > 3:  #can only interpolate with enough data points
                        cs = interp1d(depths,temps,kind='cubic')
                        xs = np.arange(depths[0],depths[len(depths)-1],1) #depths spaced by 1m
                        theta = cs(xs)
                        #HC_i = 0
                       # for k in range(len(xs)): 
                        #    HC_i += theta[k]
                        ones = np.zeros((len(xs)))
                        for i in range(len(ones)):
                            ones[i] = 1   
                        HC_grid[lt, ln] = np.dot(theta,ones) #sums over all interpolated temperatures more efficiently
                        HC_grid *= rho*Cp
    return HC_grid


def monthly_avgs(HC):
    monthlies = np.zeros((int(len(HC)/12),12))
    counter_m = 0
    counter_n = 0
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


def run(start_year, end_year, depth_from, depth_to):
    years, times, rootgrps = retrieve(1950,2018)
    
    HC = calculate_HC(rootgrps,25,31, -43, 41)
    
    months, month_avgs = monthly_avgs(HC)
    pos = str(-43)+"N "+str(41)+"E"
    
    return years, times, HC, pos, months, month_avgs

years, times, HC, pos, months, month_avgs = run(1950, 2018, 25, 31)

    
#fig, ax = plt.subplots(figsize=(6.5, 4))
#ax.plot(depths,temps0,"x")
#ax.plot(xs0,cs0(xs0))
#plt.title(str(years[0]))
#plt.grid()
#plt.show()
print(time.time()-start)

csfont = {'fontname':'Helvetica'}
hfont = {'fontname':'Helvetica'}

plt.rcParams['axes.facecolor'] = '#e3e1d8'
plt.rcParams['axes.edgecolor'] = '#e3e1d8'
fig2, ax = plt.subplots(figsize=(6.5, 4))
ax.plot(times,HC, color="#d6261a")
#plt.fill_between(times, HC-HC_Err, HC+HC_Err,
#alpha=1, edgecolor='#3F7F4C', facecolor='#7EFF99',
#linewidth=0)
plt.title("HC in 1000m-2000m depth layer from "+str(years[0])+" to "+str(years[len(years)-1])+" at \n "+pos, **hfont)
plt.xlabel("Year", **csfont)
plt.ylabel("Heat Content $(Jm^{-2})$", **csfont, labelpad= 1)
plt.grid(color = "white")
#plt.grid()
plt.show()

plt.rcParams['axes.facecolor'] = '#cccccc'
fig3, ax = plt.subplots(figsize=(6.5, 4))
ax.plot(months,month_avgs,color="#d6261a")
plt.title("Seasonal variation in HC in 1000m-2000m depth layer from "+str(years[0])+" to "+str(years[len(years)-1])+" at \n "+pos+" with monthly averages removed", **hfont)
plt.xlabel("Year", **csfont)
plt.ylabel("Heat Content Monthly Average ($Jm^{-2}$)", **csfont, labelpad = 1)
plt.grid(color="white")
plt.show()