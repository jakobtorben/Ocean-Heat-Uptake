# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:47:28 2019

@author: benro

NB: lat = lat index - 83
    lon = lon index + 1
"""
from numba import jit
import matplotlib.pyplot as plt
import numpy as np
import os
os.environ["PROJ_LIB"] = "C:/Users/jakob/Anaconda3/Library/share"; #fixr
import time
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import interp1d
import netCDF4 as nt

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
    """Input depth range in terms of indices, input lattitude and longitude (not the indices),
    calculate_errors doesn't currently work. Outputs an array of length <the number of 
    months analysed>. NB: the array named - ones - exists so that the dot product of it 
    with thetas, the temperature array, can be found: equivalent, but more efficient than 
    summing over thetas."""
    depth_from, depth_to = depth_ind(rootgrps, depth_from, depth_to)
    depths = rootgrps[0]["depth"][depth_from:depth_to]
    xs = np.arange(depths[0],depths[len(depths)-1],1) #depths spaced by 1m
    ones = np.zeros((len(xs))) #ones for the dot product with thetas (equivalent to summing over thetas)
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
    depths_i = rootgrps[0]["depth"][:]
    HC_grid_i = np.ma.zeros([len(lat), len(lon)])
    HC_grid = np.ma.masked_all_like(HC_grid_i)

    for i in range(len(rootgrps)):  #loops through every month
        for lt in range(len(lat)):
            for ln in range(len(lon)):
                temp_mask = rootgrps[i]["temperature"][0,depth_from:(depth_to+1),lt,ln].mask
                if np.ma.all(temp_mask) == False: #checks if all temp data at location are empty  
                    if np.ma.any(temp_mask) == True:   # finds data with empty temps to separate values depth<depthrange
                        pos = np.where(temp_mask == True)[0][0]
                        depths = depths_i[depth_from:(pos+1)] #finds the limited depth
                    else:
                        depths = depths_i[depth_from:(depth_to+1)] #depth is given by depthrange  
                    temps = (rootgrps[i]["temperature"][0,depth_from:(depth_from+len(depths)),lt,ln])
                    if len(temps) > 3:  #can only interpolate with enough data points
                        cs = interp1d(depths,temps, kind='cubic')
                        xs = np.arange(depths[0],depths[len(depths)-1],1) #depths spaced by 1m
                        theta = cs(xs)
                        ones = np.ones((len(xs)))
                        HC_grid[lt, ln] = np.dot(theta,ones) #sums over all interpolated temperatures more efficiently
    HC_grid *= rho*Cp
    return HC_grid


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


def run(start_year, end_year, depth_from, depth_to):
    """Basically a main; runs all functions defined above."""
    years, times, rootgrps = retrieve(1950,2018)
    
    HC = calculate_HC(rootgrps,25,31, -43, 41)
    
    months, month_avgs = monthly_avgs(HC)
    pos = str(-43)+"N "+str(41)+"E"
    
    return years, times, HC, pos, months, month_avgs

def run_global(start_year, end_year, depth_from, depth_to, animate=True):
    """Basically a main; runs all functions defined above."""
#    years, times, rootgrps = retrieve(1950,2018)
#    rootgrps_1950 = [nt.Dataset("EN.4.2.1.analyses.g10.1950\EN.4.2.1.f.analysis.g10.195001.nc", "r+", format="NETCDF4")]
#    rootgrps_1951 = [nt.Dataset("EN.4.2.1.analyses.g10.1951\EN.4.2.1.f.analysis.g10.195101.nc", "r+", format="NETCDF4")]
#    rootgrps_1952 = [nt.Dataset("EN.4.2.1.analyses.g10.1952\EN.4.2.1.f.analysis.g10.195201.nc", "r+", format="NETCDF4")]
#    rootgrps_1953 = [nt.Dataset("EN.4.2.1.analyses.g10.1953\EN.4.2.1.f.analysis.g10.195301.nc", "r+", format="NETCDF4")]
#
#
#    rootgrps_2015 = [nt.Dataset("EN.4.2.1.analyses.g10.2015\EN.4.2.1.f.analysis.g10.201501.nc", "r+", format="NETCDF4")]
#    rootgrps_2016 = [nt.Dataset("EN.4.2.1.analyses.g10.2016\EN.4.2.1.f.analysis.g10.201601.nc", "r+", format="NETCDF4")]
#    rootgrps_2017 = [nt.Dataset("EN.4.2.1.analyses.g10.2017\EN.4.2.1.f.analysis.g10.201701.nc", "r+", format="NETCDF4")]
    rootgrps_2018 = [nt.Dataset("EN.4.2.1.analyses.g10.2018\EN.4.2.1.f.analysis.g10.201801.nc", "r+", format="NETCDF4")]

#    HC_1950 = calculate_HC_global(rootgrps_1950, 0, 2000)
#    print('1950', time.time()-start)
#    HC_1951 = calculate_HC_global(rootgrps_1951, 0, 2000)
#    print('1951', time.time()-start)
#    HC_1952 = calculate_HC_global(rootgrps_1952, 0, 2000)
#    print('1952', time.time()-start)
#    HC_1953 = calculate_HC_global(rootgrps_1953, 0, 2000)
#    print('1953', time.time()-start)    
#
#    HC_2015 = calculate_HC_global(rootgrps_2015, 0, 2000)
#    print('2015', time.time()-start)
#    HC_2016 = calculate_HC_global(rootgrps_2016, 0, 2000)
#    print('2016', time.time()-start)
#    HC_2017 = calculate_HC_global(rootgrps_2017, 0, 2000)
#    print('2017', time.time()-start)
    HC_2018 = calculate_HC_global(rootgrps_2018, 0, 2000)
#    print('2018', time.time()-start)
#    HC_1950_mean = (HC_1950+HC_1951+HC_1952+HC_1953)/4
#    HC_2018_mean = (HC_2015+HC_2016+HC_2017+HC_2018)/4

#    dHC = (HC_2018_mean-HC_1950_mean)/(65*365*24*3600)
    if animate == True:
        plot(rootgrps_2018, HC_2018)
    return HC_2018

def plot(rootgrps, HC_grid):
    lon = rootgrps[0]["lon"][:]
    lat = rootgrps[0]["lat"][:]
#    clevs = np.linspace(-5, 5, 30)   #contour levels
    lon2, lat2 = np.meshgrid(lon,lat)
    fig = plt.figure(figsize=(10,7))
    m = Basemap(projection='cyl',lat_0=0,lon_0=180,resolution='l')
    x, y = m(lon2, lat2)
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawmapboundary(fill_color='white')
    m.fillcontinents('grey')
    CS = m.contourf(x, y, HC_grid, 30, cmap=plt.cm.coolwarm)
    cbar = fig.colorbar(CS, format='%.1f', shrink=0.7, orientation='horizontal')
#    cbar.ax.set_xlabel('Ocean Heat Content Change [W/m^-2]')
    plt.show()

#years, times, HC, pos, months, month_avgs = run(1950, 2018, 25, 31)
HC = run_global(1950, 2018, 0, 1000)



    
#fig, ax = plt.subplots(figsize=(6.5, 4))
#ax.plot(depths,temps0,"x")
#ax.plot(xs0,cs0(xs0))
#plt.title(str(years[0]))
#plt.grid()
#plt.show()
print(time.time()-start)
"""
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
"""