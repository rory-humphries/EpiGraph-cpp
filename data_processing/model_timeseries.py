#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 09:02:56 2020

@author: roryh
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import geopandas as gpd
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
#import imageio


files = []
for filename in os.listdir("../output/"):
    files += [filename]
    
data = pd.read_csv("../output/0.csv")

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
for filename in range(1, len(files)):
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
    
plt.plot(s, label = 'S')
plt.plot(i, label = 'I')
plt.plot(r, label = 'R')
plt.plot(d, label = 'D')
plt.plot(x, label = 'X')

plt.ylim(1, 7000000)
plt.yscale('log')
plt.legend()

#plt.annotate('local max', xy=(0.05, 0),  xycoords='axes fraction',
#            xytext=(0.05, -0.2), textcoords='axes fraction',
#            arrowprops=dict(facecolor='black', shrink=0.0005),
#           horizontalalignment='right', verticalalignment='top',
#            )
"""
fs = 8
plt.annotate('Initial Lockdown', xy=(50, 1),  xycoords='data',
            xytext=(50, 0.1), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                            headlength = 5, headwidth = 5),
            horizontalalignment='right', verticalalignment='top',
            fontsize = fs)
plt.annotate('Phase 1', xy=(117, 1),  xycoords='data',
            xytext=(117, 0.2), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                            headlength = 5, headwidth = 5),
            horizontalalignment='right', verticalalignment='top',
            fontsize = fs)
plt.annotate('Phase 2', xy=(138, 1),  xycoords='data',
            xytext=(138, 0.1), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                            headlength = 5, headwidth = 5),
            horizontalalignment='right', verticalalignment='top',
            fontsize = fs)
plt.annotate('Phase 3', xy=(159, 1),  xycoords='data',
            xytext=(159, 0.05), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                            headlength = 5, headwidth = 5),
            horizontalalignment='center', verticalalignment='top',
            fontsize = fs)
plt.annotate('Phase 4', xy=(180, 1),  xycoords='data',
            xytext=(180, 0.1), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                            headlength = 5, headwidth = 5),
            horizontalalignment='left', verticalalignment='top',
            fontsize = fs)
plt.annotate('Phase 5', xy=(201, 1),  xycoords='data',
            xytext=(201, 0.2), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                            headlength = 5, headwidth = 5),
            horizontalalignment='left', verticalalignment='top',
            fontsize = fs)
"""
plt.ylabel('No. of individuals')
plt.savefig('no_lockdown_log_scale.png', bbox_inches='tight', dpi = 300)

fig, ax = plt.subplots()
plt.plot(i, label = 'I')
plt.plot(d, label = 'D')
plt.plot(x, label = 'X')
plt.legend()

"""
fs = 8
plt.annotate('Initial Lockdown', xy=(50, 0),  xycoords='data',
            xytext=(50, -100000), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                            headlength = 5, headwidth = 5),
            horizontalalignment='right', verticalalignment='top',
            fontsize = fs)
plt.annotate('Phase 1', xy=(117, 0),  xycoords='data',
            xytext=(117, -100000), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                            headlength = 5, headwidth = 5),
            horizontalalignment='right', verticalalignment='top',
            fontsize = fs)
plt.annotate('Phase 2', xy=(138, 0),  xycoords='data',
            xytext=(138, -140000), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                            headlength = 5, headwidth = 5),
            horizontalalignment='right', verticalalignment='top',
            fontsize = fs)
plt.annotate('Phase 3', xy=(159, 0),  xycoords='data',
            xytext=(159, -180000), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                            headlength = 5, headwidth = 5),
            horizontalalignment='center', verticalalignment='top',
            fontsize = fs)
plt.annotate('Phase 4', xy=(180, 0),  xycoords='data',
            xytext=(180, -140000), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                            headlength = 5, headwidth = 5),
            horizontalalignment='left', verticalalignment='top',
            fontsize = fs)
plt.annotate('Phase 5', xy=(201, 0),  xycoords='data',
            xytext=(201, -100000), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                            headlength = 5, headwidth = 5),
            horizontalalignment='left', verticalalignment='top',
            fontsize = fs)
"""
plt.ylabel('No. of individuals')


plt.savefig('no_lockdown_lin_scale.png', bbox_inches='tight', dpi = 300)

fig, ax = plt.subplots()
plt.plot(i[:100], label = 'I')
plt.plot(d[:100], label = 'D')
plt.plot(x[:100], label = 'X')
plt.legend()


fs = 8
plt.annotate('Initial Lockdown', xy=(50, 0),  xycoords='data',
            xytext=(50, -1000), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.0005, width = 0.1, 
                            headlength = 5, headwidth = 5),
            horizontalalignment='right', verticalalignment='top',
            fontsize = fs)

plt.ylabel('No. of individuals')


plt.savefig('five_phase_model_pre_lockdown.png', bbox_inches='tight', dpi = 300)
