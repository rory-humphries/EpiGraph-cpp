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

direc = '../output/'
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

cmap = plt.cm.jet  # define the colormap
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
#cmaplist[0] = (.5, .5, .5, 1.0)

# create the new map
#cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
#    'Custom cmap', cmaplist, cmap.N)

# define the bins and normalize
bounds = [0, 5, 500, 1000, 1900, 3500, 6500]
norm = matplotlib.colors.BoundaryNorm(bounds, len(bounds))

# get max values
maxs = 0
maxi = 0
maxr = 0
maxd = 0
maxx = 0
mini = 1
for filename in range(1, 500):
    df = pd.read_csv(direc+str(filename)+".csv")
    maxs = max(maxs, max(df.S))
    maxi = max(maxi, max(df.I.to_numpy()/area))
    maxr = max(maxr, max(df.R))
    maxd = max(maxd, max(df.D))
    maxx = max(maxx, max(df.X))
    

images = []
for filename in range(1, 500):
    print(filename)
    fig, ax = plt.subplots(1, 1, figsize = (5,5))
    df = pd.read_csv(direc+str(filename)+".csv")
   
    N = df.S.to_numpy()+df.I.to_numpy()+df.X.to_numpy()+df.R.to_numpy()+df.D.to_numpy()
    ed_soa_gdf['flag'] = 1
    ed_soa_gdf['I'] = df.I.to_numpy()/area
    
    #ed_soa_gdf['flag'][ed_soa_gdf['I']<5] = 0
    #ed_soa_gdf['I'] *= 100000/N
    
    #ed_soa_gdf['I'][ed_soa_gdf['flag'] == 0] = 1
        
    
    lognorm = matplotlib.colors.LogNorm(0.0001, maxi)
    #lognorm = matplotlib.colors.Normalize(0, 100000*maxi)


    #ax2=ed_soa_gdf.plot(ax=ax, linewidth=0.0001, column = 'I', cmap = plt.get_cmap('YlOrRd', len(bounds)), norm=lognorm, legend = False)
    ax2=ed_soa_gdf.plot(ax=ax, linewidth=0.0001, column = 'I', cmap = plt.get_cmap('YlOrRd'), norm=lognorm, legend = False)
    patches = ax2.collections[0]
    cbar = plt.colorbar(patches, ax=ax2, extend='max')
    cbar.ax.set_ylabel('No. of infected per km sq', fontsize = 15)
    cbar.ax.tick_params(labelsize=15) 
    
    
    plt.title('Day ' + str(filename), fontsize = 15)
    plt.axis('off')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('tmp.png')
    images.append(imageio.imread("tmp.png"))
    plt.show()
imageio.mimsave('output.gif', images)




