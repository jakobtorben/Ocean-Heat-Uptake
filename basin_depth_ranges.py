# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 11:32:48 2019

@author: benro
"""

import numpy as np
import matplotlib.pyplot as plt

def depth_ranges(start_year, end_year, start_depth, end_depth, s):
    indexes = []
    years = []
    for i in range(end_year - start_year - int(s/6)):
        if (i + float(start_year) + s/12) %10 == 0:
            indexes.append(12*i)
            years.append(i+5+start_year+int(s/12))
    
    if end_depth<=200:
        increment = 10
    elif end_depth<=300:
        increment = 30
    elif end_depth<=1000:
        increment = 50
    depth_from = 0
    depth_to = depth_from + increment
    
    south_atlantic_depths = np.zeros((len(indexes)-1, int((end_depth-start_depth)/increment)))
    south_pacific_depths = np.zeros((len(indexes)-1, int((end_depth-start_depth)/increment)))
    north_atlantic_depths = np.zeros((len(indexes)-1, int((end_depth-start_depth)/increment)))
    north_pacific_depths = np.zeros((len(indexes)-1, int((end_depth-start_depth)/increment)))
    southern_depths = np.zeros((len(indexes)-1, int((end_depth-start_depth)/increment)))
    indian_depths = np.zeros((len(indexes)-1, int((end_depth-start_depth)/increment)))
    
    depths = np.zeros((int((end_depth-start_depth)/increment)))
    
    column = 0    
    while depth_to <= end_depth:
        times_south_pacific, total_global_south_pacific = np.loadtxt("HC_efficiency_south_pacific_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        times_south_atlantic, total_global_south_atlantic = np.loadtxt("HC_efficiency_south_atlantic_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        times_indian, total_global_indian = np.loadtxt("HC_efficiency_indian_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt", delimiter = ", ")
        times_north_pacific, total_global_north_pacific = np.loadtxt("HC_efficiency_north_pacific_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        times_north_atlantic, total_global_north_atlantic = np.loadtxt("HC_efficiency_north_atlantic_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        times_southern, total_global_southern = np.loadtxt("HC_efficiency_southern_"+str(s)+"mMA_anom_"+str(start_year)+"_"+str(end_year)+"_"+str(depth_from)+"m_"+str(depth_to)+"m.txt",delimiter = ", ")
        
        means_sa = np.zeros((len(indexes)-1))
        means_sp = np.zeros((len(indexes)-1))
        means_na = np.zeros((len(indexes)-1))
        means_np = np.zeros((len(indexes)-1))
        means_so = np.zeros((len(indexes)-1))
        means_in = np.zeros((len(indexes)-1))
        for i in range(len(indexes)-1):        
            means_sa[i] = np.mean(total_global_south_atlantic[indexes[i]:indexes[i+1]])
            means_sp[i] = np.mean(total_global_south_pacific[indexes[i]:indexes[i+1]])
            means_na[i] = np.mean(total_global_north_atlantic[indexes[i]:indexes[i+1]])
            means_np[i] = np.mean(total_global_north_pacific[indexes[i]:indexes[i+1]])
            means_so[i] = np.mean(total_global_southern[indexes[i]:indexes[i+1]])
            means_in[i] = np.mean(total_global_indian[indexes[i]:indexes[i+1]])
        
        means_sa = means_sa - means_sa[0]
        means_sp = means_sp - means_sp[0]
        means_na = means_na - means_na[0]
        means_np = means_np - means_np[0]
        means_so = means_so - means_so[0]
        means_in = means_in - means_in[0]
        
        south_atlantic_depths[:,column] = means_sa
        south_pacific_depths[:,column] = means_sp
        north_atlantic_depths[:,column] = means_na
        north_pacific_depths[:,column] = means_np
        southern_depths[:,column] = means_so
        indian_depths[:,column] = means_in
        
        depths[column] = depth_from + increment/2
#        plt.plot(years[:-1], means, "x")
#        plt.show()
        
        
        depth_from += increment
        depth_to += increment
        column += 1

    fig1, ax1 = plt.subplots(figsize=(6.5, 4))
    fig2, ax2 = plt.subplots(figsize=(6.5, 4))
    fig3, ax3 = plt.subplots(figsize=(6.5, 4))
    fig4, ax4 = plt.subplots(figsize=(6.5, 4))
    fig5, ax5 = plt.subplots(figsize=(6.5, 4))
    fig6, ax6 = plt.subplots(figsize=(6.5, 4))
    for i in range(len(indexes)-1):
        ax1.plot(depths, south_atlantic_depths[i,:], ".--", label = str(years[i]-5)+"'s")
        ax2.plot(depths, south_pacific_depths[i,:], ".--", label = str(years[i]-5)+"'s")
        ax3.plot(depths, north_atlantic_depths[i,:], ".--", label = str(years[i]-5)+"'s")
        ax4.plot(depths, north_pacific_depths[i,:], ".--", label = str(years[i]-5)+"'s")
        ax5.plot(depths, southern_depths[i,:], ".--", label = str(years[i]-5)+"'s")
        ax6.plot(depths, indian_depths[i,:], ".--", label = str(years[i]-5)+"'s")
        
    ax1.set_xlabel("Depth (m)")
    ax2.set_xlabel("Depth (m)")
    ax3.set_xlabel("Depth (m)")
    ax4.set_xlabel("Depth (m)")
    ax5.set_xlabel("Depth (m)")
    ax6.set_xlabel("Depth (m)")
    ax1.set_ylabel("Change in HC per Unit Area ($Jm^{-2}$), relative to 1960's avg")
    ax2.set_ylabel("Change in HC per Unit Area ($Jm^{-2}$), relative to 1960's avg")
    ax3.set_ylabel("Change in HC per Unit Area ($Jm^{-2}$), relative to 1960's avg")
    ax4.set_ylabel("Change in HC per Unit Area ($Jm^{-2}$), relative to 1960's avg")
    ax5.set_ylabel("Change in HC per Unit Area ($Jm^{-2}$), relative to 1960's avg")
    ax6.set_ylabel("Change in HC per Unit Area ($Jm^{-2}$), relative to 1960's avg")
    ax1.legend(title = "Decadal average \n for each depth:")
    ax2.legend(title = "Decadal average \n for each depth:")
    ax3.legend(title = "Decadal average \n for each depth:")
    ax4.legend(title = "Decadal average \n for each depth:")
    ax5.legend(title = "Decadal average \n for each depth:")
    ax6.legend(title = "Decadal average \n for each depth:")
    ax1.set_title("South Atlantic")
    ax2.set_title("South Pacific")
    ax3.set_title("North Atlantic")
    ax4.set_title("North Pacific")
    ax5.set_title("Southern Ocean")
    ax6.set_title("Indian Ocean")
    ax1.grid(linestyle = "--")
    ax2.grid(linestyle = "--")
    ax3.grid(linestyle = "--")
    ax4.grid(linestyle = "--")
    ax5.grid(linestyle = "--")
    ax6.grid(linestyle = "--")
    plt.show()
    
depth_ranges(1950, 2018, 0, 200, 12*5)        