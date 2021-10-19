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

data_path = '../data/processed/init_conditions/'


max_time = 499
#####################


files = []
for filename in os.listdir(data_path):
    files += [filename]

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

fig, ax = plt.subplots(1, 1, figsize = (8,5))


filename = files[0]
df = pd.read_csv(data_path + str(filename), index_col = False)
plt.plot(dates, df.I[:len(dates)], label = "I", c = 'tab:blue')
plt.plot(dates, df.X[:len(dates)], label = 'X', c = 'tab:green')
plt.plot(dates, df.D[:len(dates)], label = 'D', c = 'tab:red')


    
for filename in files[1:]:
    df = pd.read_csv(data_path + str(filename), index_col = False)
    plt.plot(dates, df.I[:len(dates)], c = 'tab:blue')
    plt.plot(dates, df.X[:len(dates)], c = 'tab:green')
    plt.plot(dates, df.D[:len(dates)], c = 'tab:red')

        


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

plt.legend(loc=1)

if save_figs == True:
    plt.savefig('output', bbox_inches='tight', dpi=300)

