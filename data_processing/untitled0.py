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

data_paths = ['../scenarios/compliance_99.csv',
              '../scenarios/compliance_80.csv', 
              '../scenarios/compliance_70.csv', 
              '../scenarios/compliance_40.csv']

output_name = 'compliance_comp'
fig_name = output_name + '.png'
data_name = output_name + '.csv'
title = '5 phase and 6 week lockdown'
max_time = 500
#####################


    
years = mdates.YearLocator()   # every year
months = mdates.MonthLocator()  # every month
years_fmt = mdates.DateFormatter('%Y')
months_fmt = mdates.DateFormatter('%m')

start_date = datetime.datetime(2020, 1,22)#+ datetime.timedelta(days=25)

dates = [start_date + datetime.timedelta(days=k) for k in range(max_time)]

#ax = plt.subplots(1, 1, figsize = (16/3, 9/3))
fig, ax = plt.subplots(1, 1, figsize = (15/2.5, 9/2.5))

for data in [pd.read_csv(data_path, index_col=False) for data_path in data_paths]:
    #plt.plot(dates, data.S[:max_time], label = 'S', c = 'tab:purple')
    plt.plot(dates, data.I[:max_time], label = 'I', c = 'tab:blue')
    #plt.plot(dates, data.R[:max_time], label = 'R', c = 'tab:orange')
    #plt.plot(dates, data.D[:max_time], label = 'D', c = 'tab:red')
    #plt.plot(dates, data.X[:max_time], label = 'X', c = 'tab:green')


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

plt.ylim(1, data.S[:max_time].max() + 5000000)
plt.yscale('log')

ax.legend()
plt.xlim(start_date - datetime.timedelta(days=25), start_date + datetime.timedelta(days=max_time+25))
if save_figs == True:
    plt.savefig(fig_name, bbox_inches='tight', dpi = 300)

#pd.DataFrame({'S':s,'I':i,'X':x,'R':r,'D':d}).to_csv(data_name)
