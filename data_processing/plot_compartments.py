#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 00:41:18 2020

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
save_figs = False
op_name = 'varied_mixing_x'
annotate_events = False

#####################
years = mdates.YearLocator()   # every year
months = mdates.MonthLocator()  # every month
years_fmt = mdates.DateFormatter('%Y')
months_fmt = mdates.DateFormatter('%m')
    
fig, ax = plt.subplots(1, 1)
df = pd.read_csv("../agg_output.csv")

dates = [datetime.datetime(2020, 1, 22) + datetime.timedelta(days=x) for x in range(400)]
plt.plot(dates, df.S,c='tab:purple', label = 'S', alpha = 0.5)
plt.plot(dates, df.R,c='tab:orange', label = 'R', alpha = 0.5)
plt.plot(dates, df.I, c='tab:blue', label = 'I', alpha = 0.9)
plt.plot(dates, df.X, c='tab:green', label = 'X', alpha = 0.9)
plt.plot(dates, df.D, c='tab:red', label = 'D', alpha = 0.9)

df = pd.read_csv("../22742_output.csv")


plt.plot(dates, df.S.iloc[0:400], '--', c='tab:purple')
plt.plot(dates, df.R.iloc[0:400],'--', c='tab:orange')
plt.plot(dates, df.I.iloc[0:400],'--', c='tab:blue')
plt.plot(dates, df.X.iloc[0:400],'--', c='tab:green')
plt.plot(dates, df.D.iloc[0:400],'--', c='tab:red')

ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(years_fmt)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_minor_formatter(months_fmt)
ax.tick_params(which='major', length=16)
plt.setp(ax.xaxis.get_minorticklabels(), rotation=45)
plt.legend()
#plt.savefig('compliance_diff_scale.png', bbox_inches='tight', dpi = 300)
#plt.yscale('log')

#fig, ax = plt.subplots(1, 1, figsize = (16/3, 9/3))
fig, ax = plt.subplots(1, 1)
for dirs in ['../2274', '../scenarios/dyn_lockdown_10_3']:
    files = []
    for filename in os.listdir(dir):
        files += [filename]
        
    data = pd.read_csv(dir + "/0.csv")
    
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
    for filename in range(0, len(files)):
        df = pd.read_csv(dir + '/' + str(filename)+".csv")
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
    
    dates = [datetime.datetime(2020, 1, 22) + datetime.timedelta(days=x) for x in range(len(files))]
    
    if node == 0:
        #plt.plot(dates, s, lt[node], c = 'tab:purple', alpha = at[node])
        plt.plot(dates, i, lt[node], c = 'tab:blue', alpha = at[node])
        plt.plot(dates, r, lt[node], c = 'tab:orange', alpha = at[node])
        plt.plot(dates, d, lt[node], c = 'tab:red', alpha = at[node])
        plt.plot(dates, x, lt[node], c = 'tab:green', alpha = at[node])
    if node == 1:
        #plt.plot(dates, s, lt[node], label = 'S', c = 'tab:purple',  alpha = at[node])
        plt.plot(dates, i, lt[node], label = 'I', c = 'tab:blue', alpha = at[node])
        plt.plot(dates, r, lt[node], label = 'R', c = 'tab:orange', alpha = at[node])
        plt.plot(dates, d, lt[node], label = 'D', c = 'tab:red', alpha = at[node])
        plt.plot(dates, x, lt[node], label = 'X', c = 'tab:green', alpha = at[node])
    
    node += 1
ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(years_fmt)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_minor_formatter(months_fmt)
ax.tick_params(which='major', length=16)
plt.setp(ax.xaxis.get_minorticklabels(), rotation=45)
plt.legend(loc = 0)
plt.yscale('log')
plt.ylim(1,(50e6))
    
    

    #plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
    
    #plt.xlabel('time (months)')
plt.ylabel('No. of individuals')
    
    
if save_figs == True:
    plt.savefig(op_name + '_log_scale.png', bbox_inches='tight', dpi = 300)