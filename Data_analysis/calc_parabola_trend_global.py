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


def HC_parabola_gradient(HC, times, start_year, end_year, depth_from, depth_to):
    HC_gradient = np.ma.zeros([end_year - (start_year+1) +1, 173, 360])
    HC_2dfit_error = np.ma.zeros([3, 173, 360])
    HC_parab_diff = np.ma.zeros([1, 173, 360])
    x = np.arange(start_year +1, end_year+1)
    for lt in range(173):
        print("Completed : ", lt/173*100, "%")
        for ln in range(360):
            HC_p = HC[:, lt, ln]
            times1, HC_avg_p = moving_avg(times, HC_p, 12)
            FIT, COV = np.polyfit(times1, HC_avg_p, 2 ,cov=True)
            y = x*x*FIT[0] + x*FIT[1] + FIT[2]
            y1 = (start_year+1)**2*FIT[0] + (start_year+1)*FIT[1] + FIT[2]  # +1 as moving avg removes 1 year
            y2 = end_year**2*FIT[0] + end_year*FIT[1] + FIT[2]  # -1 as moving avg removes 1 year
            HC_parab_diff[0, lt ,ln] = (y2 - y1)/((end_year - (start_year+1) +1)*3600*24*365.242) # Divide by whole time period in s
            for i in range(3):
                HC_2dfit_error[i, lt ,ln] = np.sqrt(np.diag(COV))[i]/(3600*24*365.242) #convert year to s
            for i in range(end_year - (start_year+1) +1):
                HC_gradient[i, lt ,ln] = np.gradient(y)[i]/(3600*24*365.242) #convert year to s
    return HC_gradient, HC_parab_diff, HC_2dfit_error

def validation(HC, times, start_year, end_year, depth_from, depth_to):
    HC_gradient = np.ma.zeros([end_year - (start_year+1) +1, 173, 360])
    HC_2dfit_error = np.ma.zeros([3, 173, 360])
    HC_parab_diff = np.ma.zeros([1, 173, 360])
    x = np.arange(start_year +1, end_year+1)
    for lt in range(23,143, 10):#173):
        for ln in range(335,336):#360):
            HC_p = HC[:, lt, ln]
            times1, HC_avg_p = moving_avg(times, HC_p, 12)
            FIT, COV = np.polyfit(times1, HC_avg_p, 2 ,cov=True)
            y = x*x*FIT[0] + x*FIT[1] + FIT[2]
            y1 = (start_year+1)**2*FIT[0] + (start_year+1)*FIT[1] + FIT[2]  # +1 as moving avg removes 1 year
            y2 = end_year**2*FIT[0] + end_year*FIT[1] + FIT[2]  # -1 as moving avg removes 1 year
            plt.plot(x,y)
            plt.plot(times1, HC_avg_p)
            plt.title("Location latitude: "+ str(lt-83.)+" longitude " +str(ln))
            #print("y2 = ", y2)
            #print("y1 = ", y1)
            #print("Difference = ", y2-y1)
            print("Difference in s ", lt, ln, " = ",(y2-y1)/ ((end_year - (start_year+1) +1)*3600*24*365.242))
            plt.show()
            plt.pause(5)
            plt.close()
            for i in range(3):
                HC_2dfit_error[i, lt ,ln] = np.sqrt(np.diag(COV))[i]/(3600*24*365.242) #convert year to s
            #print("Fit error: ", HC_2dfit_error[:, lt, ln])
            for i in range(end_year - (start_year+1) +1):
                HC_gradient[i, lt ,ln] = np.gradient(y)[i]/(3600*24*365.242) #convert year to s
            #print("Gradient: ", HC_gradient[:, lt, ln])
    
    
def run(start_year, end_year, depth_from, depth_to):
    years, times, rootgrps = retrieve(start_year, end_year)
    #HC = read_HC(start_year, end_year, depth_from, depth_to)
    #HC.dump("test4")
    HC = np.load("test4", allow_pickle=True)
    HC_gradient, HC_parab_diff, HC_2dfit_error = HC_parabola_gradient(HC, times, start_year, end_year, depth_from, depth_to)

    HC_gradient.dump("HC_G_gradient_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m")
    HC_2dfit_error.dump("HC_G_2dfit_error_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m")
    HC_parab_diff.dump("HC_G_parab_diff_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m")
    print('Year', start_year, 'to', end_year, 'finished')

    
def run_validation(start_year, end_year, depth_from, depth_to):
    years, times, rootgrps = retrieve(start_year, end_year)
#    HC = read_HC(start_year, end_year, depth_from, depth_to)
#    HC.dump("test4")
    HC = np.load("test4", allow_pickle=True)  # 2000-2018 0-700m
    validation(HC, times, start_year, end_year, depth_from, depth_to)


#run(2000, 2018, 0, 700)
run_validation(2000, 2018, 0, 700)

