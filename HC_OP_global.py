# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 17:47:18 2019

@author: jakob
"""

import matplotlib.pyplot as plt
import numpy as np
import os
os.environ["PROJ_LIB"] = "C:/Users/jakob/Anaconda3/Library/share"; #fixr
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import interp1d
import netCDF4 as nt
from datetime import datetime
startTime = datetime.now()


Cp = 3850 #Heat capacity of seawater
rho = 1025 #density of seawater
year = 1950 #start year
month = 1 #(start month)
years = [] #array to keep track of the years analysed (mainly for graph titles)
times = [] #for the final plot of HC against time
#rootgrps = []
rootgrps = [nt.Dataset("EN.4.2.1.analyses.g10.1950\EN.4.2.1.f.analysis.g10.195001.nc", "r+", format="NETCDF4")]
"""
for i in range(1):
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
"""

depths = rootgrps[0]["depth"][:]
depthrange = [0, 1000]    #depth range in meters
depth1 = 0                   #algorithm to find index range of depths
depth2 = 0
count = 0

for i in range(len(depths)):
    if depths[i] >= depthrange[0] and count == 0:
        depth1 = i
        count += 1
    elif depths[i] > depthrange[1] and depths[i-1] <= depthrange[1]:
        depth2 = i - 1
        
depths = rootgrps[0]["depth"][depth1:(depth2+1)]


def vol(depthtop, depthbot, lat):
    re = 6378.137*1e3 #equatorial radius earth (m)
    rtop = re-depthtop
    e = 0.08181919  #eccentricity
    dln = dlt = 1*np.pi/180
    lt = lat*np.pi/180
    dA = rtop*rtop*np.cos(lt)*(1-e*e)*dlt*dln/(1-e*e*np.sin(lt)*np.sin(lt))**2    #area
    dV = dA*(depthbot-depthtop)  #volume
    #dr = rtop - rbot
    #volume = rbot*rbot*abs(np.sin(lat*np.pi/180))*1*1*dr  #dV=r^2*sin(ϴ)*dϴ*d*ϕ*dr
    volume = rtop*rtop*abs(np.sin(lt*np.pi/180))*dln*dlt*(depthbot-depthtop)
    return dV



lon = rootgrps[0]["lon"][:]
lat = rootgrps[0]["lat"][:]
#print('lon', len(lon), lon)
#print('lat', len(lat), lat)

ln2, lt2 = np.meshgrid(lon,lat)
HC_grid_i = (np.ma.zeros(np.shape(lt2)))
HC_grid = np.ma.masked_all_like(HC_grid_i)


temps = []
for i in range(len(rootgrps)):  #loops through every month
    for lt in range(len(lat)):#range(int(lat[0]), int(lat[len(lat)-1])+1):  #latitude
        for ln in range(len(lon)):#range(int(lon[0]), int(lon[len(lon)-1])+1):  #longitutde
            if np.ma.all(rootgrps[i]["temperature"][0,depth1:(depth2+1),lt,ln].mask) == False: #checks if all temp data at location are empty  
                if np.ma.any(rootgrps[i]["temperature"][0,depth1:(depth2+1),lt,ln].mask) == True:   # finds data with empty temps to separate values depth<depthrange
                    count = 0
                    for j in rootgrps[i]["temperature"][0,depth1:(depth2+1),lt,ln].mask:
                        if j == True:
                            count += 1
                    depths = rootgrps[i]["depth"][depth1:(depth2-count+1)] #finds the limited depth
                else:
                    depths = rootgrps[i]["depth"][depth1:(depth2+1)] #depth is given by depthrange  
                temps = (rootgrps[i]["temperature"][0,depth1:(depth1+len(depths)),lt,ln])
                if len(temps) > 3:  #can only interpolate with enough data points
                    cs = interp1d(depths,temps,kind='cubic')
                    xs = np.arange(depths[0],depths[len(depths)-1],1) #depths spaced by 1m
                    theta = cs(xs)
                    HC_i = 0
                    for k in range(len(xs)): 
                        HC_i += theta[k]
                    volume = vol(depths[0], depths[len(depths)-1], lat[lt])
                    HC_grid[lt, ln] = HC_i
                    #print(lt,ln)


HC_grid *= rho*Cp 
lon2, lat2 = np.meshgrid(lon,lat)

m = Basemap(projection='cyl',lat_0=0,lon_0=180,resolution='l')
x, y = m(lon2, lat2)
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='white')
m.fillcontinents('grey')
m.contourf(x, y, HC_grid, 30, cmap=plt.cm.coolwarm)
plt.colorbar()
plt.show()

"""


position_analysed = str(rootgrps[0]["lat"][40])+"N "+str(rootgrps[0]["lon"][40])+"E" #fetches coordinates



#fig2, ax = plt.subplots(figsize=(6.5, 4))
#ax.plot(times,HC)
#plt.title("HC in 1000m-2000m depth layer from "+str(years[0])+" to "+str(years[len(years)-1])+" at "+position_analysed)
#plt.grid()
#plt.show()

"""

print('Time to run script', datetime.now() - startTime)
