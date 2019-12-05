# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 11:40:39 2019

@author: benro
"""
import numpy as np
import matplotlib.pyplot as plt

def plot_gbl_HC_trends(start_year, end_year, depth_from, depth_to, s, efficiency = False):
    if efficiency == False:
        times_south_pacific, total_global_south_pacific = np.loadtxt("HC_south_pacific_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        times_south_atlantic, total_global_south_atlantic = np.loadtxt("HC_south_atlantic_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        times_indian, total_global_indian = np.loadtxt("HC_indian_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt", delimiter = ", ")
        times_north_pacific, total_global_north_pacific = np.loadtxt("HC_north_pacific_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        times_north_atlantic, total_global_north_atlantic = np.loadtxt("HC_north_atlantic_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        times_southern, total_global_southern = np.loadtxt("HC_southern_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        y_label = "HC anomaly in J from "+str(depth_from)+"m to "+str(depth_to)+"m"
        ax.ticklabel_format(axis ='y',scilimits=(22,22))
    else:
        times_south_pacific, total_global_south_pacific = np.loadtxt("HC_efficiency_south_pacific_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        times_south_atlantic, total_global_south_atlantic = np.loadtxt("HC_efficiency_south_atlantic_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        times_indian, total_global_indian = np.loadtxt("HC_efficiency_indian_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt", delimiter = ", ")
        times_north_pacific, total_global_north_pacific = np.loadtxt("HC_efficiency_north_pacific_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        times_north_atlantic, total_global_north_atlantic = np.loadtxt("HC_efficiency_north_atlantic_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        times_southern, total_global_southern = np.loadtxt("HC_efficiency_southern_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        y_label = "HC anomaly per unit area ($Jm^{-2}$) from "+str(depth_from)+"m to "+str(depth_to)+"m"
    
    plt.rcParams['axes.facecolor'] = 'white' #cccccc
    plt.rcParams['axes.edgecolor'] = 'black'
    fig2, ax = plt.subplots(figsize=(6.5, 4))
#    plt.fill_between(times2, np.zeros((len(total_global2))), total_global2,
#    alpha=1, edgecolor='#1a4f61', facecolor='#1a4f61',
#    linewidth=0, label = "0m to 2000m (12 month moving average)")
    ax.plot(times_south_pacific, total_global_south_pacific, color = "#0072de", linewidth = "1.5", label="South Pacific")
    ax.plot(times_north_pacific, total_global_north_pacific, color = "#00386e", linewidth = "1.5", label="North Pacific")
    
#    plt.fill_between(times, np.zeros((len(total_global))), total_global,
#    alpha=1, edgecolor='#3990ad', facecolor='#3990ad',
#    linewidth=0, label = "0m to 1000m (12 month moving average)")
    ax.plot(times_south_atlantic, total_global_south_atlantic, color = "#f7260f", linewidth = "1.5", label="South Atlantic")
    ax.plot(times_north_atlantic, total_global_north_atlantic, color = "#a61b0c", linewidth = "1.5", label="North Atlantic")
    
#    plt.fill_between(times3, np.zeros((len(total_global3))), total_global3,
#    alpha=1, edgecolor='#7fb7c9', facecolor='#7fb7c9',
#    linewidth=0, label = "0m to 700m (12 month moving average)")
    ax.plot(times_indian, total_global_indian, color = "#f55511", linewidth = "1.5", label="Indian")
    ax.plot(times_southern, total_global_southern, color = "#0a8f22", linewidth = "1.5", label = "Southern")
    #ax.set_axisbelow(True)
    plt.xlabel("Year")
    plt.ylabel(y_label)

    plt.legend()
    plt.grid(linestyle = "--")
    #plt.grid(linestyle = "--", color="white")
    plt.show()

plot_gbl_HC_trends(1950,2018,600,650, 12*5, True)