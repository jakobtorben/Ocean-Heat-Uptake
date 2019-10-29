# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:01:34 2019

@author: benro
"""

import numpy as np
import netCDF4 as nt
import os
os.environ["PROJ_LIB"] = "C:/Users/benro/.conda/pkgs/proj4-5.2.0-ha925a31_1/Library/share"; #fixr
import time
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import imageio


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

def read_HC(start_year, end_year,depth_from, depth_to):
    n_years = end_year - start_year + 1
    print(n_years)
    HC = np.ma.zeros([n_years*12, 173, 360])
    year = start_year
    month = 1
    for i in range(n_years):
        for j in range(12):
            HC_j = np.loadtxt("HC_grid_"+"year_"+str(year)+"_month_"+str(month)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter=", ")
            for m in range(173):
                for n in range(360):
                    HC[(j+i*12),m,n] = HC_j[m,n]
            month += 1
        year += 1
        month = 1
                    
    return HC

def HC_glob_avg(HC):
    HC_global_avg = np.zeros((len(HC[:,0,0])))
    for i in range(len(HC[:,0,0])):
        non_zero_HC = []
        #HC_global_avg[i] = np.mean(HC[i,:,:])
        for lat in range(len(HC[i,:,0])):
            for lon in range(len(HC[i,lat,:])):
                if HC[i,lat,lon] != 0:
                    non_zero_HC.append(HC[i,lat,lon])
        non_zero_HC = np.asarray(non_zero_HC)
        HC_global_avg[i] = np.mean(non_zero_HC)
    
    return HC_global_avg

def plot_HC_maps(start_year, depth_from, depth_to, rootgrps, HC_grid):
    filenames = []
    plt.ioff()
    print(len(rootgrps))
    for i in range(len(rootgrps)):
        time = start_year+float(i)/12
        lon = rootgrps[0]["lon"][:]
        lat = rootgrps[0]["lat"][:]
        lon2, lat2 = np.meshgrid(lon,lat)
        fig = plt.figure(figsize=(10,7))
        m = Basemap(projection='cyl',lat_0=0,lon_0=180,resolution='l')
        x, y = m(lon2, lat2)
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))
        m.drawmapboundary(fill_color='white')
        m.fillcontinents('grey')
        CS = m.contourf(x, y, HC_grid[i,:,:], 30, cmap=plt.cm.coolwarm)
        cbar = fig.colorbar(CS, shrink=0.7, orientation='horizontal')
        plt.title('HC per m^2: ' +str(time))
        fig.savefig("HC_"+str(depth_from)+"m_"+str(depth_to)+"m_"+str(time)+".png")
        filenames.append("HC_"+str(depth_from)+"m_"+str(depth_to)+"m_"+str(time)+".png")
        #plt.show()
        plt.close(fig)
    
    images = []
    for filename in filenames:
        images.append(imageio.imread(filename))
    imageio.mimsave("HC_"+str(depth_from)+"m_"+str(depth_to)+"m_from_"+str(start_year)+"_to_"+str(time)+".gif",
                images, duration=0.4)
    


def run(start_year, end_year, depth_from, depth_to):
    years, times, rootgrps = retrieve(start_year, end_year)
    HC = read_HC(start_year, end_year, depth_from, depth_to)
    HC_global_avg = HC_glob_avg(HC)
    print(HC_global_avg)
    plt.plot(times, HC_global_avg)
    plt.xlabel("Year")
    plt.ylabel("Mean global HC (Jm^-2)")
    plt.grid(color="white")
    plt.show()
    #plot(start_year, depth_from, depth_to, rootgrps, HC)


run(1950, 1959, 0, 1000)
