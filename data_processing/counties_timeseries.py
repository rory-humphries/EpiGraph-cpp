#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 09:02:56 2020

@author: roryh
"""

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

data_path = '../data/output/'

output_name = 'county_lockdown_comp_100'
fig_name = output_name + '.png'
data_name = output_name + '.csv'
title = 'County Level Lockdown'

max_time = 600
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

fig, ax = plt.subplots(1, 1, figsize = (16/2.5, 9/2.5))

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
    for filename in range(0, max_time):
        df = pd.read_csv(data_path + str(filename) + ".csv")
        # s.append(sum(df.S))
        i.append(sum(df.loc[ed_soa_gdf.county == county, 'I']))
        # r.append(sum(df.R))
        # d.append(sum(df.D))
        # x.append(sum(df.X))

    # maxs = max(s);
    maxi = max(i);
    # maxr = max(r);
    # maxd = max(d);
    # maxx = max(x)

    # fig, ax = plt.subplots(1, 1, figsize = (16/3, 9/3))
    # plt.plot(dates, s, label = 'S', c = 'tab:purple')
    if county == "DUBLIN":
        l1=plt.plot(dates, i, '-', c = 'tab:blue', label=county)
    elif county == "CORK":
        l2=plt.plot(dates, i, '-', c = 'tab:green', label=county)
    elif county == "ANTRIM":
        l3=plt.plot(dates, i, '-', c = 'tab:orange', label=county)
    elif county == "DOWN":
        l4=plt.plot(dates, i, '-', c = 'darkorchid', label=county)
    else:
        plt.plot(dates, i, 'r-', alpha=0.2)
    # plt.plot(dates, r, label = 'R', c = 'tab:orange')
    # plt.plot(dates, d, label = 'D', c = 'tab:red')
    # plt.plot(dates, x, label = 'X', c = 'tab:green')
    data.append(i)


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

#plt.title(title)
plt.xlabel('Time [months]')
plt.ylabel('No. of individuals')

# plt.yscale('log')
plt.ylim(0,3500)
handles, labels = ax.get_legend_handles_labels()
ax.legend([handles[2], handles[0], handles[1], handles[3]], 
          [labels[2], labels[0], labels[1], labels[3]],
          loc=1)

if save_figs == True:
    plt.savefig(fig_name, bbox_inches='tight', dpi=300)

pd.DataFrame({c: d for c, d in zip(counties, data)}).to_csv(data_name)
