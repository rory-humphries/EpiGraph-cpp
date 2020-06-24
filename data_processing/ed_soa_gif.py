#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 10:35:36 2020

@author: roryh
"""

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

direc = '../scenarios/257_1_0/'
files = []
for filename in os.listdir(direc):
    files += [filename]
    

# align shapefile with electoral divisions by osied
# complicated as the shapefile combines some

# get the osied column for ed's as integers
ed_info = pd.read_csv('../data/raw/ED_Basic_Info.csv')
ed_osied = [int(x.split('-')[0].strip().lower()) for x in ed_info['Electoral Division'] ]
ed_info.index = ed_osied


# hold info on every Electoral Division (ed) and super output area (soa)
ed_soa_df = pd.read_csv('../data/raw/Joined_Pop_Data_CSO_NISRA.csv')
new_ind = []
for x in ed_soa_df['Electoral Division'] :
    xsplit = x.split('-')
    new_ind.append(xsplit[0].strip())
ed_soa_df.index = new_ind

# read in shapefile
ed_gdf = gpd.read_file('../data/raw/Shapefiles/electoral_divisions/')
soa_gdf = gpd.read_file('../data/raw/Shapefiles/super_output_areas')

# find which osied's are combined in the shapefile
poly_id = dict()
for x, geom in zip(ed_gdf.OSIED, ed_gdf.geometry):
    xsplit = x.split('/')
    for i in xsplit:
        if i[0] == '0':
            i=i[1:]
        poly_id[i] = geom
        
for x, geom in zip(soa_gdf.SOA_CODE, soa_gdf.geometry):
    poly_id[x] = geom
    
geoms = []
for ind in ed_soa_df.index:
    geoms.append(poly_id[ind])
    
ed_soa_gdf = gpd.GeoDataFrame(ed_soa_df)

ed_soa_gdf['geometry'] = geoms
area = ed_soa_gdf.area.to_numpy()/1000

# get max values
maxs = 0
maxi = 0
maxr = 0
maxd = 0
maxx = 0
for filename in range(1, len(files)):
    df = pd.read_csv(direc+str(filename)+".csv")
    maxs = max(maxs, max(df.S))
    maxi = max(maxi, max(df.I/area))
    maxr = max(maxr, max(df.R))
    maxd = max(maxd, max(df.D))
    maxx = max(maxx, max(df.X))
    

images = []
for filename in range(1, len(files)):
    print(filename)
    fig, ax = plt.subplots(1, 1, figsize = (10,10))
    df = pd.read_csv(direc+str(filename)+".csv")
   
    
    ed_soa_gdf['I'] = 0.1
    ed_soa_gdf['I'] = df.I.to_numpy()/area
    
    #lognorm = matplotlib.colors.LogNorm(1, max(1, df.I.max()))
    lognorm = matplotlib.colors.Normalize(1, maxi)


    ax2=ed_soa_gdf.plot(ax=ax, linewidth=0.0001, column = 'I', norm=lognorm, legend = False)
    patches = ax2.collections[0]
    cbar = plt.colorbar(patches, ax=ax2, extend='max')
    cbar.ax.set_ylabel('No. of infected', fontsize = 15)
    cbar.ax.tick_params(labelsize=15) 
    
    
    plt.title('Day ' + str(filename), fontsize = 15)
    plt.axis('off')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('tmp.png')
    images.append(imageio.imread("tmp.png"))
    plt.show()
imageio.mimsave('output.gif', images)



day = 50
fig, ax = plt.subplots(1, 1, figsize = (10,10))
df = pd.read_csv(direc+ str(day) +".csv")
   

ed_soa_gdf['I'] = 0.1
ed_soa_gdf['I'] = df.I.to_numpy()/area

lognorm = matplotlib.colors.LogNorm(0.00001, maxi)
#lognorm = matplotlib.colors.Normalize(0.00001, maxi)


ax2=ed_soa_gdf.plot(ax=ax, linewidth=0.0001, column = 'I', norm=lognorm, legend = False)
patches = ax2.collections[0]
cbar = plt.colorbar(patches, ax=ax2, extend='min')
cbar.ax.set_ylabel('No. of infected per km sq', fontsize = 15)
cbar.ax.tick_params(labelsize=15) 


plt.title('Day ' + str(day), fontsize = 15)
plt.axis('off')
plt.axis('equal')
plt.tight_layout()
plt.savefig('tmp.png')
images.append(imageio.imread("tmp.png"))
plt.show()


