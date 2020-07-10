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
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

######################
save_figs = True
op_name = 'dynamic_lock_10_3'
annotate_events = False
dir_path = "../output/"
#####################


files = []
for filename in os.listdir(dir_path):
    files += [filename]
    
num_ts = 500#len(files)

s=[]
i=[]
r=[]
d=[]
x=[]

maxs = 0
maxi = 0
maxr = 0
maxd = 0
maxx = 0
for filename in range(0, num_ts):
    df = pd.read_csv("../output/"+str(filename)+".csv")
    maxs = max(maxs, max(df.S))
    maxi = max(maxi, max(df.I))
    maxr = max(maxr, max(df.R))
    maxd = max(maxd, max(df.D))
    maxx = max(maxx, max(df.X))
    s.append(sum(df.S))
    i.append(sum(df.I))
    r.append(sum(df.R))
    d.append(sum(df.D))
    x.append(sum(df.X))
    
years = mdates.YearLocator()   # every year
months = mdates.MonthLocator()  # every month
years_fmt = mdates.DateFormatter('%Y')
months_fmt = mdates.DateFormatter('%m')

start_date = datetime.datetime(2020, 1, 22)

dates = [start_date + datetime.timedelta(days=x) for x in range(num_ts)]

#fig, ax = plt.subplots(1, 1, figsize = (16/3, 9/3))
fig, ax = plt.subplots(1, 1)
plt.plot(dates, s, label = 'S', c = 'tab:purple')
plt.plot(dates, i, label = 'I', c = 'tab:blue')
plt.plot(dates, r, label = 'R', c = 'tab:orange')
plt.plot(dates, d, label = 'D', c = 'tab:red')
plt.plot(dates, x, label = 'X', c = 'tab:green')

event_durations = [50, 67, 21, 21, 21, 21, 21, 400]
event_dates = [start_date]
for x in event_durations:
    event_dates.append(event_dates[-1] + datetime.timedelta(days=x))
    
shade_event_dates = event_dates[1:-1]
num_shade_colors = len(shade_event_dates)

colors = plt.cm.viridis(np.linspace(0,1,num_shade_colors))

for i in range(len(shade_event_dates)-1):
    d1 = shade_event_dates[i]
    d2 = shade_event_dates[i+1]
    plt.axvspan(d1, d2, color = colors[i], alpha=0.2)

ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(years_fmt)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_minor_formatter(months_fmt)
ax.tick_params(which='major', length=16)
plt.setp(ax.xaxis.get_minorticklabels(), rotation=45)
#plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)

#plt.xlabel('time (months)')
plt.ylabel('No. of individuals')

plt.ylim(1, 9000000)
plt.yscale('log')
plt.legend()

if save_figs == True:
    plt.savefig(op_name + '_log_scale.png', bbox_inches='tight', dpi = 300)

"""

if annotate_events == True:
    fs = 8
    diff = 0.4
    p1 = 0.05
    p2 = 0.15
    plt.annotate('Initial Lockdown', xy=(datetime.datetime(2020, 3, 12), 1),  xycoords='data',
                xytext=(datetime.datetime(2020, 3, 12), p1), textcoords='data',
                arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                                headlength = 5, headwidth = 5),
                horizontalalignment='right', verticalalignment='top',
                fontsize = fs)
    plt.annotate('Phase 1', xy=(datetime.datetime(2020, 5, 18), 1),  xycoords='data',
                xytext=(datetime.datetime(2020, 5, 18), p2), textcoords='data',
                arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                                headlength = 5, headwidth = 5),
                horizontalalignment='right', verticalalignment='top',
                fontsize = fs)
    plt.annotate('Phase 2', xy=(datetime.datetime(2020, 6, 8), 1),  xycoords='data',
                xytext=(datetime.datetime(2020, 6, 8), p2*diff), textcoords='data',
                arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                                headlength = 5, headwidth = 5),
                horizontalalignment='right', verticalalignment='top',
                fontsize = fs)
    plt.annotate('Phase 3', xy=(datetime.datetime(2020, 6, 29), 1),  xycoords='data',
                xytext=(datetime.datetime(2020, 6, 29), p2*diff*diff), textcoords='data',
                arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                                headlength = 5, headwidth = 5),
                horizontalalignment='center', verticalalignment='top',
                fontsize = fs)
    plt.annotate('Phase 4', xy=(datetime.datetime(2020, 7, 20), 1),  xycoords='data',
                xytext=(datetime.datetime(2020, 7, 20),  p2*diff), textcoords='data',
                arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                                headlength = 5, headwidth = 5),
                horizontalalignment='left', verticalalignment='top',
                fontsize = fs)
    plt.annotate('Phase 5', xy=(datetime.datetime(2020, 8, 10), 1),  xycoords='data',
                xytext=(datetime.datetime(2020, 8, 10), p2), textcoords='data',
                arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                                headlength = 5, headwidth = 5),
                horizontalalignment='left', verticalalignment='top',
                fontsize = fs)

if save_figs == True:
    plt.savefig(op_name + '_log_scale.png', bbox_inches='tight', dpi = 300)


#fig, ax = plt.subplots(1, 1, figsize = (16/3, 9/3))
fig, ax = plt.subplots(1, 1)
plt.plot(dates, i, label = 'I', c = 'tab:blue')
plt.plot(dates, d, label = 'D', c = 'tab:red')
plt.plot(dates, x, label = 'X', c = 'tab:green')
plt.legend()

if annotate_events == True:
    fs = 8
    diff = -25000
    p1 = -120000
    p2 = -100000
    fs = 8
    plt.annotate('Initial Lockdown', xy=(datetime.datetime(2020, 3, 12), 0),  xycoords='data',
                xytext=(datetime.datetime(2020, 3, 12), p1), textcoords='data',
                arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                                headlength = 5, headwidth = 5),
                horizontalalignment='right', verticalalignment='top',
                fontsize = fs)
    plt.annotate('Phase 1', xy=(datetime.datetime(2020, 5, 18), 0),  xycoords='data',
                xytext=(datetime.datetime(2020, 5, 18), p2), textcoords='data',
                arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                                headlength = 5, headwidth = 5),
                horizontalalignment='right', verticalalignment='top',
                fontsize = fs)
    plt.annotate('Phase 2', xy=(datetime.datetime(2020, 6, 8), 0),  xycoords='data',
                xytext=(datetime.datetime(2020, 6, 8), p2+diff), textcoords='data',
                arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                                headlength = 5, headwidth = 5),
                horizontalalignment='right', verticalalignment='top',
                fontsize = fs)
    plt.annotate('Phase 3', xy=(datetime.datetime(2020, 6, 29), 0),  xycoords='data',
                xytext=(datetime.datetime(2020, 6, 29), p2 + 2*diff), textcoords='data',
                arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                                headlength = 5, headwidth = 5),
                horizontalalignment='center', verticalalignment='top',
                fontsize = fs)
    plt.annotate('Phase 4', xy=(datetime.datetime(2020, 7, 20), 0),  xycoords='data',
                xytext=(datetime.datetime(2020, 7, 20), p2 + diff), textcoords='data',
                arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                                headlength = 5, headwidth = 5),
                horizontalalignment='left', verticalalignment='top',
                fontsize = fs)
    plt.annotate('Phase 5', xy=(datetime.datetime(2020, 8, 10), 0),  xycoords='data',
                xytext=(datetime.datetime(2020, 8, 10), p2), textcoords='data',
                arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                                headlength = 5, headwidth = 5),
                horizontalalignment='left', verticalalignment='top',
                fontsize = fs)

ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(years_fmt)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_minor_formatter(months_fmt)
ax.tick_params(which='major', length=16)
plt.setp(ax.xaxis.get_minorticklabels(), rotation=45)
#plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)

#plt.xlabel('time (months)')
plt.ylabel('No. of individuals')

if save_figs == True:
    plt.savefig(op_name + '_lin_scale.png', bbox_inches='tight', dpi = 300)

"""
