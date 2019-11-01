# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 11:40:39 2019

@author: benro
"""
import numpy as np
import matplotlib.pyplot as plt

def plot_gbl_HC_trends(start_year, end_year, depth_from, depth_to, moving_average = True, anomaly_plot = True):
    if moving_average == True:
        if anomaly_plot == True:
            times, total_global = np.loadtxt("HC_G_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
            y_label = "Total global HC anomaly "+str(depth_from)+"m-"+str(depth_to)+"m in J (moving average)"
        else:
            times, total_global = np.loadtxt("HC_G_MA_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
            y_label = "Total global HC "+str(depth_from)+"m-"+str(depth_to)+"m in J (moving average)"
    else:
        if anomaly_plot == True:
            times, total_global = np.loadtxt("HC_G_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
            y_label = "Total global HC anomaly "+str(depth_from)+"m-"+str(depth_to)+"m in J"
        else:
            times, total_global = np.loadtxt("HC_G_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
            y_label = "Total global HC "+str(depth_from)+"m-"+str(depth_to)+"m in J"
    
    plt.rcParams['axes.facecolor'] = 'white' #cccccc
    plt.rcParams['axes.edgecolor'] = 'black'
    fig2, ax = plt.subplots(figsize=(6.5, 4))
    plt.fill_between(times, np.zeros((len(total_global))), total_global,
    alpha=0.8, edgecolor='#d99b00', facecolor='#d99b00',
    linewidth=0)
    ax.plot(times, total_global, color = "#b86e00")
    #ax.set_axisbelow(True)
    plt.xlabel("Year")
    plt.ylabel(y_label)
    #plt.grid(linestyle = "--", color="white")
    plt.show()

plot_gbl_HC_trends(1950,2018,0,1000)