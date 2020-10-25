# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:01:34 2019

@author: benro
"""

import numpy as np
import netCDF4 as nt
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
    HC = np.ma.zeros([n_years*12, 173, 360])
    year = start_year
    month = 1
    for i in range(n_years):
        if i%5==0:
            print("Read "+str(i)+" of "+str(n_years)+" years")
        for j in range(12):
            HC_j = np.loadtxt("HC_grid_"+"year_"+str(year)+"_month_"+str(month)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter=", ")
            for m in range(173):
                for n in range(360):
                    HC[(j+i*12),m,n] = HC_j[m,n]
            month += 1
        year += 1
        month = 1
                    
    return HC

def read_HC_uncertainty(start_year, end_year,depth_from, depth_to):
    n_years = end_year - start_year + 1
    HC = np.ma.zeros([n_years*12, 173, 360])
    year = start_year
    month = 1
    for i in range(n_years):
        if i%5==0:
            print("Read "+str(i)+" of "+str(n_years)+" years")
        for j in range(12):
            HC_j = np.loadtxt("HC_grid_uncertainty_"+"year_"+str(year)+"_month_"+str(month)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter=", ")
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
        if i%60 == 0:
            print("Calculated "+str(i/12)+" of "+str(len(HC[:,0,0])/12)+" years")
        for lat in range(len(HC[0,:,0])):
            dA = area(depth_from,lat-83)
            for lon in range(len(HC[0,0,:])):
                ttl_global[i] += HC[i,lat,lon]*dA
    return ttl_global
                
    

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
    


def run(start_year, end_year, depth_from, depth_to):
    years, times, rootgrps = retrieve(start_year, end_year)
    HC = read_HC(start_year, end_year, depth_from, depth_to)
#    HC_global_avg = HC_glob_avg(HC)
#    times1, HC_global_avg = moving_avg(times, HC_global_avg, 12)
    plt.rcParams['axes.facecolor'] = '#e3e1d8'
    plt.rcParams['axes.edgecolor'] = '#e3e1d8'
#    fig1, ax = plt.subplots(figsize=(6.5, 4))
#    ax.plot(times1, HC_global_avg)
#    plt.xlabel("Year")
#    plt.ylabel("Moving Average of (Unweighted) Mean global HC 0m-1000m (Jm^-2)")
#    plt.grid(color="white")
#    plt.show()
    
    total_global = ttl_global_HC(depth_from, HC)
    np.savetxt("HC_G_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times,total_global),delimiter=", ")
    np.savetxt("HC_G_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times,std_baseline(total_global)),delimiter=", ")
    times2, total_global = moving_avg(times, total_global, 12)
    np.savetxt("HC_G_MA_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2, total_global),delimiter=", ")
    total_global = std_baseline(total_global)
    np.savetxt("HC_G_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global),delimiter=", ")

def run_uncertainty(start_year, end_year, depth_from, depth_to):
    years, times, rootgrps = retrieve(start_year, end_year)
    HC = read_HC(start_year, end_year, depth_from, depth_to)
    HC_uncrt = read_HC_uncertainty(start_year, end_year, depth_from, depth_to)

    
    total_global_plus = ttl_global_HC(depth_from, HC + HC_uncrt)
    np.savetxt("HC_G_uncertainty_plus_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times,total_global_plus),delimiter=", ")
    np.savetxt("HC_G_anom_uncertainty_plus_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times,std_baseline(total_global_plus)),delimiter=", ")
    
    total_global_minus = ttl_global_HC(depth_from, HC - HC_uncrt)
    np.savetxt("HC_G_uncertainty_minus_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times,total_global_minus),delimiter=", ")
    np.savetxt("HC_G_anom_uncertainty_minus_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times,std_baseline(total_global_minus)),delimiter=", ")
    
    times2, total_global_plus = moving_avg(times, total_global_plus, 12)
    np.savetxt("HC_G_MA_uncertainty_plus_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2, total_global_plus),delimiter=", ")
    total_global_plus = std_baseline(total_global_plus)
    np.savetxt("HC_G_MA_anom_uncertainty_plus_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_plus),delimiter=", ")

    times2, total_global_minus = moving_avg(times, total_global_minus, 12)
    np.savetxt("HC_G_MA_uncertainty_minus_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2, total_global_minus),delimiter=", ")
    total_global_minus = std_baseline(total_global_minus)
    np.savetxt("HC_G_MA_anom_uncertainty_minus_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",(times2,total_global_minus),delimiter=", ")

    plt.rcParams['axes.facecolor'] = '#e3e1d8'
    plt.rcParams['axes.edgecolor'] = '#e3e1d8'
    fig1, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(times2, total_global_plus)
    ax.plot(times2, total_global_minus)
    plt.xlabel("Year")
    plt.ylabel("Moving Average of (Unweighted) Mean global HC 0m-1000m (Jm^-2)")
    plt.grid(color="white")
    plt.show()

run(1950, 2018, 0, 2000)
run(1950, 2018, 0, 1000)
run(1950, 2018, 0, 700)
run(1950, 2018, 0, 100)
run(1950, 2018, 0, 300)
#run_uncertainty(1950, 2018, 0, 2000)