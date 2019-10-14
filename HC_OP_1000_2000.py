# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 18:25:20 2019

@author: benro
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


Cp = 3850 #Heat capacity of seawater
rho = 1025 #density of seawater
year = 1950 #start year
month = 1 #(start month)
years = [] #array to keep track of the years analysed (mainly for graph titles)
times = [] #for the final plot of HC against time
HC = []

rootgrps = []

for i in range(69):
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


depths = rootgrps[0]["depth"][25:31] #selects depth range from roughly 950 - 2100m
xs = np.arange(970,2120,1) #depths spaced by 1m

for i in range(len(rootgrps)):
    temps = rootgrps[i]["temperature"][0,25:31,40,40]
    cs = interp1d(depths,temps,kind='cubic')
    theta = cs(xs)
    HC_i = 0
    for j in range(len(xs)):
        HC_i += theta[j]
    HC.append(HC_i)
HC  = np.asarray(HC)
HC *= rho*Cp


position_analysed = str(rootgrps[0]["lat"][40])+"N "+str(rootgrps[0]["lon"][40])+"E" #fetches coordinates
    
#fig, ax = plt.subplots(figsize=(6.5, 4))
#ax.plot(depths,temps0,"x")
#ax.plot(xs0,cs0(xs0))
#plt.title(str(years[0]))
#plt.grid()
#plt.show()

fig2, ax = plt.subplots(figsize=(6.5, 4))
ax.plot(times,HC)
plt.title("HC in 1000m-2000m depth layer from "+str(years[0])+" to "+str(years[len(years)-1])+" at "+position_analysed)
plt.grid()
plt.show()