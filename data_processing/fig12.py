import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cbook as cbook

import datetime

import os
import geopandas as gpd
from shapely import wkt

######################
save_figs = True
max_time = 600
#####################

# read in shapefile
ed_soa_gdf = pd.read_csv('../data/processed/ed_soa_data_frame.csv')
ed_soa_gdf['geometry'] = ed_soa_gdf['geometry'].apply(wkt.loads)
ed_soa_gdf = gpd.GeoDataFrame(ed_soa_gdf, crs='epsg:4326')
ed_soa_gdf = ed_soa_gdf.to_crs('epsg:29902')

counties = ed_soa_gdf.county.unique()

start_date = datetime.datetime(2020, 1, 22)

dates = [start_date + datetime.timedelta(days=x) for x in range(max_time)]
years = mdates.YearLocator()  # every year
months = mdates.MonthLocator()  # every month
years_fmt = mdates.DateFormatter('%Y')
months_fmt = mdates.DateFormatter('%m')

fig, ax = plt.subplots(1, 1, figsize = (16/2.5, 9/2.5))

dirs = ["../scenarios/comp_70/", 
        "../scenarios/comp_80/",
        "../scenarios/comp_90/",
        "../scenarios/comp_99/",
        "../scenarios/comp_100/"]
labels = ['70%', '80%', '90%', '99%', '100%']

for data_path, label in zip(dirs, labels):
    data = []
    for county in counties:
        print(county)
    
        s = [];
        i = [];
        r = [];
        d = [];
        x = []
    
        maxs = 0;
        maxi = 0;
        maxr = 0;
        maxd = 0;
        maxx = 0
        if county != "DUBLIN":
            continue
        
        for filename in range(0, max_time):
            df = pd.read_csv(data_path + str(filename) + ".csv")
            i.append(sum(df.loc[ed_soa_gdf.county == county, 'I']))

        maxi = max(i);
        plt.plot(dates, i, '-', label=label)

        data.append(i)

plt.legend()
plt.hlines(1400, dates[0], dates[-1], color='k', linestyle = '--')

event_durations = [50, 15, 52, 21, 21, 21, 21, 21, 400]
event_dates = [start_date]
for x in event_durations:
    event_dates.append(event_dates[-1] + datetime.timedelta(days=x))

shade_event_dates = event_dates[1:-1]
num_shade_colors = len(shade_event_dates)

colors = plt.cm.viridis(np.linspace(0, 1, num_shade_colors))

for i in range(len(shade_event_dates) - 1):
    d1 = shade_event_dates[i]
    d2 = shade_event_dates[i + 1]
    plt.axvspan(d1, d2, color=colors[i], alpha=0.2)

ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(years_fmt)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_minor_formatter(months_fmt)
ax.tick_params(which='major', length=16)
plt.setp(ax.xaxis.get_minorticklabels(), rotation=45)


plt.xlabel('Time [months]')
plt.ylabel('No. of individuals')
plt.ylim(0,3500)


if save_figs == True:
    plt.savefig('output', bbox_inches='tight', dpi=600)

