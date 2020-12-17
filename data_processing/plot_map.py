#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 23:59:41 2020

@author: roryh
"""



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import datetime

import os
import geopandas as gpd
from shapely import wkt

#######################################
data_path = '../scenarios/12_mar_01_dec_hist_proj_after/'
output_name = '12_mar_01_dec_hist_proj_after_map1'

cmap = 'viridis'
time = 50+15+52+21+21#+21+21+21+45+47+45+47
cb_axis_title = r'No. of infected per km$^2$'

agg_by_county = False
per_area = True
#######################################

files = []
for filename in os.listdir(data_path):
    files += [filename]
    
start_date = datetime.datetime(2020, 1, 22)
date = start_date + datetime.timedelta(days=time)

# read in shapefile
ed_soa_gdf = pd.read_csv('../data/processed/ed_soa_data_frame.csv')
ed_soa_gdf['geometry'] = ed_soa_gdf['geometry'].apply(wkt.loads)
ed_soa_gdf = gpd.GeoDataFrame(ed_soa_gdf, crs='epsg:4326')
ed_soa_gdf = ed_soa_gdf.to_crs('epsg:29902')

ed_soa_gdf['area'] = ed_soa_gdf.area/1000

if agg_by_county:
    ed_soa_gdf['area'] = ed_soa_gdf.groupby('county')['area'].transform(sum)



fig, ax = plt.subplots(1, 1, figsize = (7,7))
df = pd.read_csv(data_path+str(time)+".csv")
   
N = df.N
ed_soa_gdf['flag'] = 1
ed_soa_gdf['I'] = df.I.to_numpy()

if agg_by_county:
    ed_soa_gdf['I'] = ed_soa_gdf.groupby('county')['I'].transform(sum)
    
if per_area:
    ed_soa_gdf['I'] /= ed_soa_gdf['area']


#max_i = max(ed_soa_gdf.I)
lognorm = matplotlib.colors.LogNorm(0.000001, max_i)
#lognorm = matplotlib.colors.Normalize(0, 200)


#ax2=ed_soa_gdf.plot(ax=ax, linewidth=0.0001, column = 'I', cmap = plt.get_cmap('YlOrRd', len(bounds)), norm=lognorm, legend = False)
ax2=ed_soa_gdf.plot(ax=ax, linewidth=0.0001, column = 'I', cmap = plt.get_cmap(cmap), norm=lognorm, legend = False)
patches = ax2.collections[0]
cbar = plt.colorbar(patches, ax=ax2, extend='max')
cbar.ax.set_ylabel(cb_axis_title, fontsize = 15)
cbar.ax.tick_params(labelsize=15) 


plt.title(date.strftime("%Y-%m-%d"), fontsize=20)
plt.axis('off')
plt.axis('equal')
plt.tight_layout()

plt.savefig(output_name)





