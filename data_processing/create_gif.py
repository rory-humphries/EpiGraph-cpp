
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
import datetime
from tqdm import tqdm

import os
import geopandas as gpd
import imageio
from shapely import wkt

#######################################
data_path = '../output/'
gif_name = 'output.gif'

cmap = 'RdYlGn_r'
max_time = 500
cb_axis_title = 'No. of infected per county'

agg_by_county = True
per_area = False
output_gif = True
#######################################

files = []
for filename in os.listdir(data_path):
    files += [filename]
    
start_date = datetime.datetime(2020, 1, 22)
dates = [start_date + datetime.timedelta(days=x) for x in range(max_time)]

# read in shapefile
ed_soa_gdf = pd.read_csv('../data/processed/ed_soa_data_frame.csv')
ed_soa_gdf['geometry'] = ed_soa_gdf['geometry'].apply(wkt.loads)
ed_soa_gdf = gpd.GeoDataFrame(ed_soa_gdf, crs='epsg:4326')
ed_soa_gdf = ed_soa_gdf.to_crs('epsg:29902')

ed_soa_gdf['area'] = ed_soa_gdf.area/1000

if agg_by_county:
    ed_soa_gdf['area'] = ed_soa_gdf.groupby('county')['area'].transform(sum)


# get max values
maxs = 0;maxi = 0;maxr = 0;maxd = 0;maxx = 0;mini = 1
for filename in range(1, max(max_time, len(files))):
    df = pd.read_csv(data_path+str(filename)+".csv")
    maxs = max(maxs, max(df.S))
    maxi = max(maxi, max(df.I))
    maxr = max(maxr, max(df.R))
    maxd = max(maxd, max(df.D))
    maxx = max(maxx, max(df.X))
    

images = []
for filename in tqdm(range(1, min(max_time, len(files)))):
    fig, ax = plt.subplots(1, 1, figsize = (7,7))
    df = pd.read_csv(data_path+str(filename)+".csv")
   
    N = df.N
    ed_soa_gdf['flag'] = 1
    ed_soa_gdf['I'] = df.I.to_numpy()
    
    if agg_by_county:
        ed_soa_gdf['I'] = ed_soa_gdf.groupby('county')['I'].transform(sum)
        
    if per_area:
        ed_soa_gdf['I'] /= ed_soa_gdf['area']
    
    
    
    #lognorm = matplotlib.colors.LogNorm(1, maxi)
    lognorm = matplotlib.colors.Normalize(0, 200)


    #ax2=ed_soa_gdf.plot(ax=ax, linewidth=0.0001, column = 'I', cmap = plt.get_cmap('YlOrRd', len(bounds)), norm=lognorm, legend = False)
    ax2=ed_soa_gdf.plot(ax=ax, linewidth=0.0001, column = 'I', cmap = plt.get_cmap(cmap), norm=lognorm, legend = False)
    patches = ax2.collections[0]
    cbar = plt.colorbar(patches, ax=ax2, extend='max')
    cbar.ax.set_ylabel(cb_axis_title, fontsize = 15)
    cbar.ax.tick_params(labelsize=15) 
    
    
    plt.title(dates[filename].strftime("%Y-%m-%d"), fontsize=20)
    plt.axis('off')
    plt.axis('equal')
    plt.tight_layout()
    if output_gif:
        plt.savefig('tmp.png')
        images.append(imageio.imread("tmp.png"))
    #plt.show()
    plt.close(fig)
    
if output_gif:
    imageio.mimsave('output.gif', images[0:499])




