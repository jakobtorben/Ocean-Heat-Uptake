import numpy as np
from netCDF4 import Dataset
import os
os.environ["PROJ_LIB"] = "C:/Users/jakob/Anaconda3/Library/share"; #fixr
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import imageio


d_1960_01 = Dataset('EN.4.2.1.f.analysis.g10.196001.nc')
depth = d_1960_01.variables['depth'][:]
lat = d_1960_01.variables['lat'][:]
lon = d_1960_01.variables['lon'][:]
time = d_1960_01.variables['time'][:]
temperature = d_1960_01.variables['temperature'][:]
temperature_uncertainty = d_1960_01.variables['temperature_uncertainty'][:]

depthrange = [2000, 5000]    #depth range in meters
depth1 = 0                   #algorithm to find index range of depths
depth2 = 0
count = 0
for i in range(0, len(depth)):
    if depth[i] >= depthrange[0] and count == 0:
        depth1 = i
        count += 1
    elif depth[i] > depthrange[1] and depth[i-1] <= depthrange[1]:
        depth2 = i - 1


lon2, lat2 = np.meshgrid(lon,lat)
filenames = []
plt.ioff()
clevs = np.linspace(np.min(temperature[0,depth1,:,:]),
                    np.max(temperature[0, depth2,:,:]), 20)   #contour levels

for i in np.linspace(depth1, depth2, depth2-depth1, dtype=int):
    fig = plt.figure(figsize=(10,7))
    m = Basemap(projection='cyl',lat_0=0,lon_0=180,resolution='l')
    x, y = m(lon2, lat2)
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawmapboundary(fill_color='white')
    m.fillcontinents('grey')
    CS = m.contourf(x,y, temperature[0,i,:,:], clevs, cmap=plt.cm.coolwarm)
    cbar = fig.colorbar(CS, format='%.1f', shrink=0.7)
    cbar.ax.set_ylabel('Potential Temperature [K]')    #chagne to °C
    print('temperature = ', np.mean(temperature[0,i,:,:]))
    print('depth = ',  depth[i])
    plt.title('Temperature [K] at Depth = %1.0f m' %depth[i])
    fig.savefig('Depth%1.0f.png' %depth[i])
    filenames.append('Depth%1.0f.png' %depth[i])
    plt.close(fig)
    

images = []
for filename in filenames:
    images.append(imageio.imread(filename))
imageio.mimsave('Depth%1.0fm_to_%1.0fm.gif' %(depth[depth1] ,depth[depth2]),
                images, duration=0.4)




