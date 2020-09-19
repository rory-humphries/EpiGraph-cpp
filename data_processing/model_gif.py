#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 09:02:56 2020

@author: roryh
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import geopandas as gpd
import imageio
files = []
for filename in os.listdir("../output/"):
    files += [filename]
    
# get max values
maxs = 0
maxi = 0
maxr = 0
maxd = 0
maxx = 0
for filename in range(1, len(files)):
    df = pd.read_csv("../output/"+str(filename)+".csv")
    maxs = max(maxs, max(df.S))
    maxi = max(maxi, max(df.I))
    maxr = max(maxr, max(df.R))
    maxd = max(maxd, max(df.D))
    maxx = max(maxx, max(df.X))
    
# align shapefile with electoral divisions by osied
# complicated as the shapefile combines some

# get the osied column for ed's as integers
ed_info = pd.read_csv('../data/raw/ED_Basic_Info.csv')
ed_osied = [int(x.split('-')[0].strip().lower()) for x in ed_info['Electoral Division'] ]
ed_info.index = ed_osied

# read in shapefile
gdf = gpd.read_file('../data/raw/Shapefiles/electoral_divisions/')

# find which osied's are combined in the shapefile
gdf_osied_combined = []
for x in gdf.OSIED:
    xsplit = x.split('/')
    if len(xsplit) > 1:
        gdf_osied_combined.append([int(y) for y in xsplit])
new_index = {x[1]:x[0] for x in gdf_osied_combined}

gdf_osied = []
for x in gdf.OSIED:
    xsplit = x.split('/')
    gdf_osied.append(int(xsplit[0]))

gdf.index = gdf_osied

images = []
for filename in range(1, len(files)):
    print(filename)
    fig, ax = plt.subplots(1, 1, figsize = (7,7))
    df = pd.read_csv("../output/"+str(filename)+".csv")
    df.index = ed_info.index
    df = df.rename(index = new_index)
    
    gdf['I'] = 0.1
    gdf['I'].loc[df.index] = df.I.to_numpy()
    
    #lognorm = matplotlib.colors.LogNorm(1, max(1, df.I.max()))
    lognorm = matplotlib.colors.Normalize(1, maxi)


    ax2=gdf.plot(ax=ax, linewidth=0.0001, column = 'I', norm=lognorm, legend = False)
    patches = ax2.collections[0]
    cbar = plt.colorbar(patches, ax=ax2, extend='max')
    cbar.ax.set_ylabel('No. of infected', fontsize = 15)
    cbar.ax.tick_params(labelsize=15) 
    
    
    plt.title('Day ' + str(filename), fontsize = 15)
    plt.axis('off')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(filename+'tmp.png')
    images.append(imageio.imread("tmp.png"))
    plt.show()
#imageio.mimsave('output.gif', images)
"""
fig, ax = plt.subplots(1, 1, figsize = (10,10))
df = pd.read_csv("../output/"+str(50)+".csv")
df.index = ed_info.index
df = df.rename(index = new_index)

gdf['I'] = 0.1
gdf['I'].loc[df.index] = df.I.to_numpy()

#lognorm = matplotlib.colors.LogNorm(1, max(1, df.I.max()))
lognorm = matplotlib.colors.Normalize(1, maxi)


ax2=gdf.plot(ax=ax, linewidth=0.0001, column = 'I', norm=lognorm, legend = False)
patches = ax2.collections[0]
cbar = plt.colorbar(patches, ax=ax2, extend='max')
cbar.ax.set_ylabel('No. of infected', fontsize = 15)
cbar.ax.tick_params(labelsize=15) 


plt.title('Day ' + str(filename), fontsize = 15)
plt.axis('off')
plt.axis('equal')
plt.tight_layout()
plt.savefig(filename+'tmp.png')
plt.show()
"""