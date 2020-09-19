#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 10:35:50 2020

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
op_name = '1000_lockdown'
annotate_events = False
#####################

num_ts = 600#len(files)

fig, ax = plt.subplots(1, 1)
lables = ['40%','80%']
count = 0
for direc in ['../../scenarios/compliance_0_4', '../../scenarios/compliance_0_8']:
    files = []
    for filename in os.listdir(direc):
        files += [filename]
        
    
    s=[]
    i=[]
    r=[]
    d=[]
    x=[]

    for filename in range(0, num_ts):
        df = pd.read_csv(direc + '/' + str(filename)+".csv")
        s.append(sum(df.S))
        i.append(sum(df.I))
        r.append(sum(df.R))
        d.append(sum(df.D))
        x.append(sum(df.X))
        
    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator()  # every month
    years_fmt = mdates.DateFormatter('%Y')
    months_fmt = mdates.DateFormatter('%m')
    
    dates = [datetime.datetime(2020, 1, 22) + datetime.timedelta(days=x) for x in range(num_ts)]
    #plt.plot(dates, s, label = 'S', c = 'tab:purple')
    #plt.plot(dates, i, label = 'I', c = 'tab:blue')
   # plt.plot(dates, i, label = lables[count])
    #plt.plot(dates, d, label = 'D', c = 'tab:red')
    #plt.plot(dates, x, label = 'X', c = 'tab:green')

    
    if count == 0:
        #plt.plot(dates, s, lt[node], c = 'tab:purple', alpha = at[node])
        plt.plot(dates, i, '--', c = 'tab:blue', alpha = 0.4)
        plt.plot(dates, r, '--', c = 'tab:orange', alpha = 0.4)
        plt.plot(dates, d, '--', c = 'tab:red', alpha = 0.4)
        plt.plot(dates, x, '--', c = 'tab:green', alpha = 0.4)
    if count == 1:
        #plt.plot(dates, s, lt[node], label = 'S', c = 'tab:purple',  alpha = at[node])
        plt.plot(dates, i, label = 'I', c = 'tab:blue')
        plt.plot(dates, r, label = 'R', c = 'tab:orange')
        plt.plot(dates, d, label = 'D', c = 'tab:red')
        plt.plot(dates, x, label = 'X', c = 'tab:green')
    
    count+=1
ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(years_fmt)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_minor_formatter(months_fmt)
ax.tick_params(which='major', length=16)
plt.setp(ax.xaxis.get_minorticklabels(), rotation=45)
plt.legend(loc = 0)
plt.yscale('log')
plt.ylim(1,900000)
    
    

    #plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
    
    #plt.xlabel('time (months)')
plt.ylabel('No. of individuals')
plt.xlabel('Time [months]')

    

    
if save_figs == True:
    plt.savefig('varied_compliance.png', bbox_inches='tight', dpi = 300)

""" 
node = 0
lt = ['-', '-', '-', '-', '-']
at = [1, 1, 1, 1, 1]
#fig, ax = plt.subplots(1, 1, figsize = (16/3, 9/3))
fig, ax = plt.subplots(1, 1)
for dir in ['../scenarios/r1', '../scenarios/r2', '../scenarios/r3',  
            '../scenarios/r4',  '../scenarios/r5']:
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
    
    if node > 0:
        plt.plot(dates, s, lt[node], c = 'tab:purple', alpha = at[node])
        plt.plot(dates, i, lt[node], c = 'tab:blue', alpha = at[node])
        plt.plot(dates, r, lt[node], c = 'tab:orange', alpha = at[node])
        plt.plot(dates, d, lt[node], c = 'tab:red', alpha = at[node])
        plt.plot(dates, x, lt[node], c = 'tab:green', alpha = at[node])
    else:
        plt.plot(dates, s, lt[node], label = 'S', c = 'tab:purple',  alpha = at[node])
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