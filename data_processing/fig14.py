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

import datetime
import os

######################
data_path = '../scenarios/comp_80/'
event_durations = [49, 15, 52, 21, 21, 21, 21, 21, 45, 47, 45, 47, 400]

data_path = '../scenarios/comp_80/'
max_time = 600
#####################


files = []
for filename in os.listdir(data_path):
    files += [filename]
    


s=[];i=[];r=[];d=[];x=[]

maxs = 0;maxi = 0;maxr = 0;maxd = 0;maxx = 0
for filename in range(0, max_time):
    df = pd.read_csv(data_path+str(filename)+".csv")
    s.append(sum(df.S))
    i.append(sum(df.I))
    r.append(sum(df.R))
    d.append(sum(df.D))
    x.append(sum(df.X))
    
maxs = max(s)
maxi = max(i)
maxr = max(r)
maxd = max(d)
maxx = max(x)
    
years = mdates.YearLocator()   # every year
months = mdates.MonthLocator()  # every month
years_fmt = mdates.DateFormatter('%Y')
months_fmt = mdates.DateFormatter('%m')

start_date = datetime.datetime(2020, 1,22)#+ datetime.timedelta(days=25)

dates = [start_date + datetime.timedelta(days=k) for k in range(max_time)]

#ax = plt.subplots(1, 1, figsize = (16/3, 9/3))
fig, ax = plt.subplots(1, 1, figsize = (16/2.5, 9/2.5))
plt.plot(dates, s, label = 'S', c = 'tab:purple')
plt.plot(dates, np.array(i), label = 'I', c = 'tab:blue')
plt.plot(dates, r, label = 'R', c = 'tab:orange')
plt.plot(dates, d, label = 'D', c = 'tab:red')
plt.plot(dates, x, label = 'X', c = 'tab:green')



event_dates = [start_date]
for k in event_durations:
    event_dates.append(event_dates[-1] + datetime.timedelta(days=k))
    
shade_event_dates = event_dates[1:-1]
num_shade_colors = len(shade_event_dates)

colors = plt.cm.viridis(np.linspace(0,1,num_shade_colors))

for k in range(len(shade_event_dates)-1):
    if k in [7, 9]:
        continue
    d1 = shade_event_dates[k]
    d2 = shade_event_dates[k+1]
    plt.axvspan(d1, d2, color = colors[k], alpha=0.2)

ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(years_fmt)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_minor_formatter(months_fmt)
ax.tick_params(which='major', length=16)
plt.setp(ax.xaxis.get_minorticklabels(), rotation=45)

plt.xlabel('Time [months]')
plt.ylabel('No. of individuals')

plt.ylim(1, maxs + 1000000)
plt.yscale('log')
plt.legend(loc=4)

if save_figs == True:
    plt.savefig('output', bbox_inches='tight', dpi = 300)

