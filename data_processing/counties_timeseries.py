#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:47:57 2020

@author: roryh
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
import os

######################
save_figs = True
op_name = 'counties'
dir_path = "../output/"
#####################

files = []
for filename in os.listdir(dir_path):
    files += [filename]
    
num_ts = 500#len(files)
years = mdates.YearLocator()   # every year
months = mdates.MonthLocator()  # every month
years_fmt = mdates.DateFormatter('%Y')
months_fmt = mdates.DateFormatter('%m')
start_date = datetime.datetime(2020, 1, 22)
dates = [start_date + datetime.timedelta(days=x) for x in range(num_ts)]


df = pd.read_csv("../output/0.csv")
inds = df.sort_values('N', ascending = False).index[1:6]

# hold info on every Electoral Division (ed) and super output area (soa)
ed_soa_df = pd.read_csv('../data/raw/Joined_Pop_Data_CSO_NISRA.csv')

# maps each ed to and soa to it's index in ed_soa_df
ed_soa_id_map = {ed:x for ed, x in zip(ed_soa_df['Electoral Division'], ed_soa_df.index)}

soa_counties = pd.read_csv('../data/raw/Geographic Data (statistical geographies).csv')
soa_counties.index = soa_counties['SOA Code']

ed_counties = pd.read_csv('../data/raw/ED_Basic_Info.csv')
ed_counties.index = ed_counties['Electoral Division'].str.split(expand=True)[0]

county_set = []
county_set.extend(np.unique(ed_counties['COUNTY'].to_numpy()))
county_dict_rev = {y:x for x,y in zip(county_set, range(len(county_set)))}

plot_data = {i:[] for i in inds}

maxi = 0
for filename in range(0, num_ts):
    df = pd.read_csv("../output/"+str(filename)+".csv")
    for i in inds:
        plot_data[i].append(df.iloc[i].I)
        maxi = max(maxi, plot_data[i][-1])
     
fig, ax = plt.subplots(1, 1)

colors = plt.cm.viridis(np.linspace(0,1,len(inds)))
c = 0
for i in inds:
    plt.plot(dates, plot_data[i], label = county_dict_rev[i], c = colors[c])
    c+=1


ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(years_fmt)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_minor_formatter(months_fmt)
ax.tick_params(which='major', length=16)
plt.setp(ax.xaxis.get_minorticklabels(), rotation=45)
plt.legend()
plt.yscale('log')
plt.ylim(1, maxi)

if save_figs == True:
    plt.savefig(op_name + '.png', bbox_inches='tight', dpi = 300)

#plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)