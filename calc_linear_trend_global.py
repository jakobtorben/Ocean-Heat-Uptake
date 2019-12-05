# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 11:54:27 2019

@author: jakob
"""

import numpy as np
import netCDF4 as nt
import matplotlib.pyplot as plt

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

def moving_avg(times, HC,s):
    m = np.zeros((len(HC)-2*s))
    times = times[s:len(m)+s]
    for i in range(len(m)):
        n = i+s
        m[i] = np.mean(HC[n-s:n+s])
        
    return times, m

def calc_HC_slope(HC, times, start_year, end_year, depth_from, depth_to):
    HC_slope = np.ma.zeros([1, 173, 360])
    HC_slope_error = np.ma.zeros([1, 173, 360]) 
    for lt in range(173):
        for ln in range(360):
            HC_p = HC[:, lt, ln]
            times1, HC_avg_p = moving_avg(times, HC_p, 12)
            FIT, COV = np.polyfit(times1, HC_avg_p, 1 ,cov=True)
            HC_slope[0, lt ,ln] = FIT[0]/(3600*24*365.242) #convert year to s
            HC_slope_error[0, lt ,ln] = np.sqrt(np.diag(COV))[0]/(3600*24*365.242) #convert year to s
    return HC_slope, HC_slope_error


def calc_HC_slope_anim1(HC, times, s, start_year, end_year, depth_from, depth_to):
    n = (end_year-start_year+1)  # s = time period (months) to find gradient for
    print(n)
    HC_slope = np.ma.zeros([n-2, 173, 360])
    HC_slope_error = np.ma.zeros([n-2, 173, 360]) 

    for lt in range(173):
        for ln in range(360):
            HC_p = HC[:, lt, ln]
            times1, HC_avg_p = moving_avg(times, HC_p, 12)
            for i in np.arange(0, n*s-2*s, s):  # n-2 as moving avg removes 2 years
                print(i, lt)
                HC_avg_p2 = HC_avg_p[i : i +s +1]
                times2 = times1[i : i +s +1]
                FIT, COV = np.polyfit(times2, HC_avg_p2, 1 ,cov=True)
                HC_slope[i//s, lt ,ln] = FIT[0]/(3600*24*365.242) #convert year to s
                HC_slope_error[i//s, lt ,ln] = np.sqrt(np.diag(COV))[0]/(3600*24*365.242) #convert year to s
    return HC_slope, HC_slope_error

def calc_HC_slope_anim_validation(HC, times, s, start_year, end_year, depth_from, depth_to):
    n = (end_year-start_year+1)*s  # s = time period (months) to find gradient for
    print(n)
    HC_slope = np.ma.zeros([n-2, 173, 360])
    HC_slope_error = np.ma.zeros([n-2, 173, 360]) 

    for lt in range(80,81):
        for ln in range(335,360):
            HC_p = HC[:, lt, ln]
            times1, HC_avg_p = moving_avg(times, HC_p, 12)
            for n in np.arange(0, n-2*s, s):  # n-2 as moving avg removes 2 years
                print(n, lt)
                HC_avg_p2 = HC_avg_p[n : n +s +1]
                times2 = times1[n : n +s +1]
                FIT, COV = np.polyfit(times2, HC_avg_p2, 1 ,cov=True)
                fit=np.poly1d(FIT)
                x = np.linspace(times2[0], times2[-1], 100)
                plt.plot(x, fit(x))
                plt.plot(times2, HC_avg_p2)
                plt.show()
                plt.pause(5)
                HC_slope[n, lt ,ln] = FIT[0]/(3600*24*365.242) #convert year to s
                HC_slope_error[n, lt ,ln] = np.sqrt(np.diag(COV))[0]/(3600*24*365.242) #convert year to s

def validation(HC, times, s, start_year, end_year, depth_from, depth_to):
    HC_slope = np.ma.zeros([1, 173, 360])
    HC_slope_error = np.ma.zeros([1, 173, 360]) 
    x = np.arange(start_year +1, end_year+1)
    for lt in range(23,143, 10):#173):
        for ln in range(335,336):#360):
            HC_p = HC[:, lt, ln]
            times1, HC_avg_p = moving_avg(times, HC_p, 12)
            FIT, COV = np.polyfit(times1, HC_avg_p, 1 ,cov=True)
            fit = np.poly1d(FIT)
            print(lt, ln, FIT[0]/(3600*24*365.242))
            HC_slope[0, lt ,ln] = FIT[0]/(3600*24*365.242) #convert year to s
            HC_slope_error[0, lt ,ln] = np.sqrt(np.diag(COV))[0]/(3600*24*365.242) #convert year to s
            plt.plot(x, fit(x))
            plt.plot(times1, HC_avg_p)
            plt.title("Location latitude: "+ str(lt-83.)+" longitude " +str(ln))
            plt.show()
            plt.pause(5)
            plt.close()
    
    
def run(start_year, end_year, depth_from, depth_to):
    years, times, rootgrps = retrieve(start_year, end_year)
    HC = read_HC(start_year, end_year, depth_from, depth_to)
    HC_slope, HC_slope_error = calc_HC_slope(HC, times, start_year, end_year, depth_from, depth_to)

    HC_slope.dump("HC_G_slope_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m")
    HC_slope_error.dump("HC_G_slope_error_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m")
    print('Year', start_year, 'to', end_year, 'finished')
    
def run_animation(s, start_year, end_year, depth_from, depth_to):
    years, times, rootgrps = retrieve(start_year, end_year)
    #HC = read_HC(start_year, end_year, depth_from, depth_to)
    #HC.dump("test1950to2018_animation")
    HC = np.load("test1950to2018_animation", allow_pickle=True)  # 2000-2018 0-700m
    #calc_HC_slope_anim_validation(HC, times, start_year, end_year, depth_from, depth_to)
    HC_slope, HC_slope_error = calc_HC_slope_anim1(HC, times, s, start_year, end_year, depth_from, depth_to)

    HC_slope.dump("HC_G_slope_animation_1year_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m")
    HC_slope_error.dump("HC_G_slope_error_animation_1year_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m")
    print('Year', start_year, 'to', end_year, 'finished')
    
def run_validation(s, start_year, end_year, depth_from, depth_to):
    years, times, rootgrps = retrieve(start_year, end_year)
    #HC = read_HC(start_year, end_year, depth_from, depth_to)
    #HC.dump("test2000to2018")
    HC = np.load("test2000to2018", allow_pickle=True)  # 2000-2018 0-700m
    validation(HC, times, s, start_year, end_year, depth_from, depth_to)

run(2000, 2018, 0, 700)
#run_validation(12, 2000, 2018, 0, 700)
#run_animation(12, 1950, 2018, 0, 2000)

