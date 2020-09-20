#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 21:29:29 2020

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
op_name = '5_phase_r0.csv'
dir_path = "../output/"
#####################

R0 = np.loadtxt("../R0.csv", delimiter = ',')


files = []
for filename in os.listdir(dir_path):
    files += [filename]
    
num_ts = 1000  # len(files)

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
#l1=plt.plot(dates, s, label = 'S', c = 'tab:purple')
l2=plt.plot(dates, i, label = 'I', c = 'tab:blue')
l3=plt.plot(dates, r, label = 'R', c = 'tab:orange')
l4=plt.plot(dates, d, label = 'D', c = 'tab:red')
l5=plt.plot(dates, x, label = 'X', c = 'tab:green')
#plt.plot(dates[0], [0], 'k--', label = r"$R_0$")

event_durations = [50, 67, 21, 21, 21, 21, 21, 400]
event_dates = [start_date]
for x in event_durations:
    event_dates.append(event_dates[-1] + datetime.timedelta(days=x))
    
shade_event_dates = event_dates[1:-1]
num_shade_colors = len(shade_event_dates)

colors = plt.cm.viridis(np.linspace(0,1,num_shade_colors))

for i in range(len(shade_event_dates) - 1):
    d1 = shade_event_dates[i]
    d2 = shade_event_dates[i + 1]
    plt.axvspan(d1, d2, color=colors[i], alpha=0.2)

ax.set_xlabel('Time [months]')
ax.set_ylabel('No. of individuals')

# plt.ylim(1, 9000000)
# plt.yscale('log')
plt.legend()

ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

ax2.xaxis.set_major_locator(years)
ax2.xaxis.set_major_formatter(years_fmt)
ax2.xaxis.set_minor_locator(months)
ax2.xaxis.set_minor_formatter(months_fmt)
ax2.tick_params(which='major', length=16)
plt.setp(ax.xaxis.get_minorticklabels(), rotation=45)
#plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
ax2.set_ylabel(r'$R_0$')  # we already handled the x-label with ax1
l6 = ax2.plot(dates, R0[:num_ts], 'k--', label=r'$R_0$', alpha=0.2)
# added these three lines
lns = l2 + l3 + l4 + l5 + l6
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=1)
#ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()