# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:01:34 2019

@author: benro

NB: lat = lat index - 83
    lon = lon index + 1
"""

import numpy as np
import os
os.environ["PROJ_LIB"] = "C:/Users/benro/.conda/pkgs/proj4-5.2.0-ha925a31_1/Library/share"; #fixr
import time
import matplotlib.pyplot as plt
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
            month += 1
        year+=1
        month = 1 #reset to the first month of the year
    return years, times, rootgrps

def read_HC(start_year, end_year,depth_from, depth_to):
    n_years = end_year - start_year + 1
    HC = np.ma.zeros([n_years*12, 173, 360])
    year = start_year
    month = 1
    for i in range(n_years):
        if i%5==0:
            print("Read "+str(i)+" of "+str(n_years)+" years")
        for j in range(12):
            HC_j = np.loadtxt("HC_grid_"+"year_"+str(year)+"_month_"+str(month)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter=", ")
#            for m in range(173):
#                for n in range(360):
#                    HC[(j+i*12),m,n] = HC_j[m,n]
            HC[j+i*12,:,:] = HC_j
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


def basins(HC, NS = False):
    atlantic_HC = np.zeros(HC.shape)
    pacific_HC = np.zeros(HC.shape)
    indian_HC = np.zeros(HC.shape)
    southern_HC = np.zeros(HC.shape)
    
    south_atlantic_HC = np.zeros(HC.shape)
    south_pacific_HC = np.zeros(HC.shape)
    north_atlantic_HC = np.zeros(HC.shape)
    north_pacific_HC = np.zeros(HC.shape)
    
    for lon in range(len(HC[0,0,:])):
        for lat in range(len(HC[0,:,0])):
            if lat-83<-30:
                southern_HC[:,lat,lon] = HC[:,lat,lon]
            else:
                if lon+1>260 or lon+1<20:
                    if lat-83>46:
                        atlantic_HC[:,lat,lon] = HC[:,lat,lon]
                    elif lat-83<30 and lat-83>=20:
                        atlantic_HC[:,lat,lon] = HC[:,lat,lon]
                    elif lat-83<20 and lat-83>8:
                        if lon+1>=260-2*(lat-83-20) or lon+1<20:
                            atlantic_HC[:,lat,lon] = HC[:,lat,lon]
                    elif lat-83<8:
                        if lon+1>292 or lon+1<20:
                            atlantic_HC[:,lat,lon] = HC[:,lat,lon]
                    elif lon+1>292 and lon+1<354:
                        atlantic_HC[:,lat,lon] = HC[:,lat,lon]
                if lon+1>=20 and lon+1<292:
                    if lon+1>20 and lon+1<117 and lat-83<30:
                        indian_HC[:,lat,lon] = HC[:,lat,lon]
                    elif lon+1>=117 and lon+1<292:
                        if lat-83>=20 and lon+1<260:
                            pacific_HC[:,lat,lon] = HC[:,lat,lon]
                        elif lat-83<20 and lat-83>=8 and lon+1<260-2*(lat-83-20):
                            pacific_HC[:,lat,lon] = HC[:,lat,lon]
                        elif lat-83<8:
                            pacific_HC[:,lat,lon] = HC[:,lat,lon]
    
    if NS == True:
        for lon in range(len(HC[0,0,:])):
            for lat in range(len(HC[0,:,0])):
                if lat-83<0:
                    south_atlantic_HC[:,lat,lon] = atlantic_HC[:,lat,lon]
                    south_pacific_HC[:,lat,lon]  = pacific_HC[:,lat,lon]
                else:
                    north_atlantic_HC[:,lat,lon] = atlantic_HC[:,lat,lon]
                    north_pacific_HC[:,lat,lon] = pacific_HC[:,lat,lon]
        return south_pacific_HC, south_atlantic_HC, north_pacific_HC, north_atlantic_HC, indian_HC, southern_HC
    else:
        return pacific_HC, atlantic_HC, indian_HC, southern_HC

    
def area(depthtop, lat):
    re = 6378.137*1e3 #equatorial radius earth (m)
    rtop = re-depthtop
    e = 0.08181919  #eccentricity
    dln = dlt = 1*np.pi/180
    lt = lat*np.pi/180
    dA = rtop*rtop*np.cos(lt)*(1-e*e)*dlt*dln/(1-e*e*np.sin(lt)*np.sin(lt))**2    #area
    return dA


def ttl_global_HC(depth_from, HC):
    ttl_global = np.zeros((len(HC[:,0,0])))
    for i in range(len(HC[:,0,0])):
        if i%180 == 0:
            print("Calculated "+str(i/12)+" of "+str(len(HC[:,0,0])/12)+" years")
#        for lat in range(len(HC[0,:,0])):
#            dA = area(depth_from,lat-83)
#            for lon in range(len(HC[0,0,:])):
#                ttl_global[i] += HC[i,lat,lon]*dA
        ttl_global[i] = np.sum(HC[i,:,:])
    return ttl_global
                

def find_ocean_area(depth_from, HC):
    ocean_area = 0
    for lat in range(len(HC[0,:,0])):
        for lon in range(len(HC[0,0,:])):
            if HC[0,lat,lon] > 1e-12:
                ocean_area += area(depth_from, lat-83)
    print(ocean_area)
    return ocean_area


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
    


def run(start_year, end_year, depth_from, depth_to, s):
    years, times, rootgrps = retrieve(start_year, end_year)
    HC = read_HC(start_year, end_year, depth_from, depth_to)
#    plt.rcParams['axes.facecolor'] = '#e3e1d8'
#    plt.rcParams['axes.edgecolor'] = '#e3e1d8'
    
    #pacific_HC, atlantic_HC, indian_HC = basins(HC)
    south_pacific_HC, south_atlantic_HC, north_pacific_HC, north_atlantic_HC, indian_HC, southern_HC = basins(HC, NS=True)
    total_global_south_pacific = ttl_global_HC(depth_from, south_pacific_HC)
    total_global_south_atlantic = ttl_global_HC(depth_from, south_atlantic_HC)
    total_global_indian = ttl_global_HC(depth_from, indian_HC)
    total_global_north_pacific = ttl_global_HC(depth_from, north_pacific_HC)
    total_global_north_atlantic = ttl_global_HC(depth_from, north_atlantic_HC)
    total_global_southern = ttl_global_HC(depth_from, southern_HC)
    
    #total_global = ttl_global_HC(depth_from, HC)
    
    #np.savetxt("HC_G_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times,total_global),delimiter=", ")

    times2, total_global_south_pacific = moving_avg(times, total_global_south_pacific, s)
    times2, total_global_south_atlantic = moving_avg(times, total_global_south_atlantic, s)
    times2, total_global_indian = moving_avg(times, total_global_indian, s)
    times2, total_global_north_pacific = moving_avg(times, total_global_north_pacific, s)
    times2, total_global_north_atlantic = moving_avg(times, total_global_north_atlantic, s)
    times2, total_global_southern = moving_avg(times, total_global_southern, s)
    
    #times2, total_global = moving_avg(times, total_global_southern, s)
    
    #np.savetxt("HC_G_MA_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2, total_global),delimiter=", ")
    #total_global = std_baseline(total_global)

    total_global_south_pacific = total_global_south_pacific - np.mean(total_global_south_pacific[:629])
    total_global_south_atlantic = total_global_south_atlantic - np.mean(total_global_south_atlantic[:629])
    total_global_indian = total_global_indian - np.mean(total_global_indian[:629])
    total_global_north_pacific = total_global_north_pacific - np.mean(total_global_north_pacific[:629])
    total_global_north_atlantic = total_global_north_atlantic - np.mean(total_global_north_atlantic[:629])
    total_global_southern = total_global_southern - np.mean(total_global_southern[:629])
    
    #total_global = total_global - np.mean(total_global_southern[:629])
#    np.savetxt("HC_south_pacific_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_south_pacific),delimiter=", ")
#    np.savetxt("HC_south_atlantic_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_south_atlantic),delimiter=", ")
#    np.savetxt("HC_indian_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_indian),delimiter=", ")
#    np.savetxt("HC_north_pacific_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_north_pacific),delimiter=", ")
#    np.savetxt("HC_north_atlantic_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_north_atlantic),delimiter=", ")
#    np.savetxt("HC_southern_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_southern),delimiter=", ")
        
        
    
    total_global_south_pacific = total_global_south_pacific/find_ocean_area(depth_from, south_pacific_HC)
    total_global_south_atlantic = total_global_south_atlantic/find_ocean_area(depth_from, south_atlantic_HC)
    total_global_indian = total_global_indian/find_ocean_area(depth_from, indian_HC)
    total_global_north_pacific = total_global_north_pacific/find_ocean_area(depth_from, north_pacific_HC)
    total_global_north_atlantic = total_global_north_atlantic/find_ocean_area(depth_from, north_atlantic_HC)
    total_global_southern = total_global_southern/find_ocean_area(depth_from, southern_HC)
    
    #total_global = total_global/find_ocean_area(depth_from, HC)
    
    np.savetxt("HC_efficiency_south_pacific_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_south_pacific),delimiter=", ")
    np.savetxt("HC_efficiency_south_atlantic_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_south_atlantic),delimiter=", ")
    np.savetxt("HC_efficiency_indian_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_indian),delimiter=", ")
    np.savetxt("HC_efficiency_north_pacific_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_north_pacific),delimiter=", ")
    np.savetxt("HC_efficiency_north_atlantic_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_north_atlantic),delimiter=", ")
    np.savetxt("HC_efficiency_southern_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_southern),delimiter=", ")
    
    #np.savetxt("HC_efficiency_global_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global),delimiter=", ")




run(1950, 2018, 0, 50, 12*5)
run(1950, 2018, 50, 100, 12*5)
run(1950, 2018, 100, 150, 12*5)
run(1950, 2018, 150, 200, 12*5)
run(1950, 2018, 200, 250, 12*5)
run(1950, 2018, 250, 300, 12*5)
run(1950, 2018, 300, 350, 12*5)
run(1950, 2018, 350, 400, 12*5)
run(1950, 2018, 400, 450, 12*5)
run(1950, 2018, 450, 500, 12*5)
run(1950, 2018, 500, 550, 12*5)
run(1950, 2018, 550, 600, 12*5)
run(1950, 2018, 600, 650, 12*5)
run(1950, 2018, 650, 700, 12*5)
run(1950, 2018, 700, 750, 12*5)
run(1950, 2018, 750, 800, 12*5)
run(1950, 2018, 800, 850, 12*5)
run(1950, 2018, 850, 900, 12*5)
run(1950, 2018, 900, 950, 12*5)
run(1950, 2018, 950, 1000, 12*5)
