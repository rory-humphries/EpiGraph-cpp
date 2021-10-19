#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 15:03:09 2020

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
import matplotlib.dates as mdates

import datetime

import os

######################
save_figs = True

data_path = '../data/output/MinImpactLDZones.csv'

output_name = 'doomsday_r0_1'
fig_name = output_name + '.png'
data_name = output_name + '.csv'
title = '5 phase and 6 week lockdown'
max_time = 600
#####################

data = pd.read_csv(data_path, index_col=False)
    
years = mdates.YearLocator()   # every year
months = mdates.MonthLocator()  # every month
years_fmt = mdates.DateFormatter('%Y')
months_fmt = mdates.DateFormatter('%m')

start_date = datetime.datetime(2020, 1,22)#+ datetime.timedelta(days=25)

dates = [start_date + datetime.timedelta(days=k) for k in range(max_time)]

#ax = plt.subplots(1, 1, figsize = (16/3, 9/3))
fig, ax = plt.subplots(1, 1, figsize = (15/2.5, 9/2.5))
#l1=plt.plot(dates, data.S[:max_time], label = 'S', c = 'tab:purple')
l2=plt.plot(dates, data[:max_time], label = 'I', c = 'tab:blue')
#l3=plt.plot(dates, data.R[:max_time], label = 'R', c = 'tab:orange')
#l4=plt.plot(dates, data.D[:max_time], label = 'D', c = 'tab:red')
#l5=plt.plot(dates, data.X[:max_time], label = 'X', c = 'tab:green')


event_durations = [50, 15, 52, 21, 21, 21, 21, 21, 400]

event_dates = [start_date]
for k in event_durations:
    event_dates.append(event_dates[-1] + datetime.timedelta(days=k))
    
shade_event_dates = event_dates[1:-1]
num_shade_colors = len(shade_event_dates)

colors = plt.cm.viridis(np.linspace(0,1,num_shade_colors))

for k in range(len(shade_event_dates)-1):
    d1 = shade_event_dates[k]
    d2 = shade_event_dates[k+1]
    plt.axvspan(d1, d2, color = colors[k], alpha=0.2)

ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(years_fmt)
ax.xaxis.set_minor_locator(months)
ax.xaxis.set_minor_formatter(months_fmt)
ax.tick_params(which='major', length=16)
plt.setp(ax.xaxis.get_minorticklabels(), rotation=45)
#plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)

#plt.title(title)
plt.xlabel('Time [months]')
plt.ylabel('No. of individuals')
plt.savefig('countytravels.png')

plt.ylim(1, data.S[:max_time].max() + 5000000)
plt.yscale('log')

ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

ax2.xaxis.set_major_locator(years)
ax2.xaxis.set_major_formatter(years_fmt)
ax2.xaxis.set_minor_locator(months)
ax2.xaxis.set_minor_formatter(months_fmt)
ax2.tick_params(which='major', length=16)
plt.setp(ax.xaxis.get_minorticklabels(), rotation=45)

tmp = [1.82, 1.28, 0.73, 0.73, 0.91, 1.09, 1.19, 1.28, 1.37]

tx = []
ty = []
for i, x in enumerate(tmp):
    tx.extend([x,x])
    ty.extend([event_dates[i], event_dates[i+1]])
    #ax2.hlines(x, event_dates[i], event_dates[i+1], linestyle = 'dotted')
l7 = ax2.plot(ty, tx, 'r-.', label = r"$R_0$")

ax2.set_ylabel(r'$R_0$')  # we already handled the x-label with ax1
l6 = ax2.plot(dates, data.R0[:max_time], 'k--', label=r'$R_0^*$', alpha=0.8)



lns = l2+l6+l7#l1+l2+l3+l4+l5+l6
labs = [l.get_label() for l in lns]
ax2.legend(lns, labs, loc=1)
plt.xlim(start_date - datetime.timedelta(days=25), start_date + datetime.timedelta(days=max_time+25))
if save_figs == True:
    plt.savefig(fig_name, bbox_inches='tight', dpi = 300)

#pd.DataFrame({'S':s,'I':i,'X':x,'R':r,'D':d}).to_csv(data_name)
