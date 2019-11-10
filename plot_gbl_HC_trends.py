# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 11:40:39 2019

@author: benro
"""
import numpy as np
import matplotlib.pyplot as plt

def plot_gbl_HC_trends(start_year, end_year, depth_from, depth_to, multi = False, moving_average = True, anomaly_plot = True):
    if moving_average == True:
        if anomaly_plot == True:
            times, total_global = np.loadtxt("HC_G_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
            y_label = "Global HC anomaly "+str(depth_from)+"m-"+str(depth_to)+"m in J (moving average)"
            if multi == True:
                times2, total_global2 = np.loadtxt("HC_G_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(2000)+"m.txt",delimiter = ", ")
                times3, total_global3 = np.loadtxt("HC_G_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(700)+"m.txt", delimiter = ", ")
                times4, total_global4 = np.loadtxt("HC_G_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(300)+"m.txt", delimiter = ", ")
                times5, total_global5 = np.loadtxt("HC_G_MA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(100)+"m.txt", delimiter = ", ")
                y_label = "Global HC anomaly in J"
        else:
            multi = False
            times, total_global = np.loadtxt("HC_G_MA_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(2000)+"m.txt",delimiter = ", ")
            y_label = "Global HC "+str(depth_from)+"m-"+str(depth_to)+"m in J (moving average)"
    else:
            times, total_global = np.loadtxt("HC_G_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
            y_label = "Global HC "+str(depth_from)+"m-"+str(depth_to)+"m in J"
    
    plt.rcParams['axes.facecolor'] = 'white' #cccccc
    plt.rcParams['axes.edgecolor'] = 'black'
    fig2, ax = plt.subplots(figsize=(6.5, 4))
    if multi == True:
        plt.fill_between(times2, np.zeros((len(total_global2))), total_global2,
        alpha=1, edgecolor='#1a4f61', facecolor='#1a4f61',
        linewidth=0, label = "0m to 2000m")
        ax.plot(times2, total_global2, color = "black", linewidth = "1.2")
        
    plt.fill_between(times, np.zeros((len(total_global))), total_global,
    alpha=1, edgecolor='#3990ad', facecolor='#3990ad',
    linewidth=0, label = "0m to 1000m")
    ax.plot(times, total_global, color = "black", linewidth = "0.7")
    
    if multi == True:
        plt.fill_between(times3, np.zeros((len(total_global3))), total_global3,
        alpha=1, edgecolor='#7fb7c9', facecolor='#7fb7c9',
        linewidth=0, label = "0m to 700m")
        ax.plot(times3, total_global3, color = "black", linewidth = "1.2")
        
        plt.fill_between(times4, np.zeros((len(total_global4))), total_global4,
        alpha=1, edgecolor='#9fc5d1', facecolor='#9fc5d1',
        linewidth=0, label = "0m to 300m")
        ax.plot(times4, total_global4, color = "black", linewidth = "1.2")
        
        plt.fill_between(times5, np.zeros((len(total_global5))), total_global5,
        alpha=1, edgecolor='#c1dde6', facecolor='#c1dde6',
        linewidth=0, label = "0m to 100m")
        ax.plot(times5, total_global5, color = "black", linewidth = "1.2")
    #ax.set_axisbelow(True)
    plt.xlabel("Year")
    plt.ylabel(y_label)
    ax.ticklabel_format(axis ='y',scilimits=(22,22))
    plt.legend(title="(12 month moving average)")
    #plt.grid(linestyle = "--", color="white")
    plt.show()

plot_gbl_HC_trends(1950,2018,0,1000, multi = True, moving_average = True, anomaly_plot = True)