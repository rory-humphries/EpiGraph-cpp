#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 11:21:40 2020

@author: roryh
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
######################
save_figs = True

output_name = 'distance_distribution'
fig_name = output_name + '.png'
data_name = output_name + '.csv'
#####################

path_list = ['../data/processed/1000.txt',
             '../data/processed/2.txt',
             '../data/processed/5.txt',
             '../data/processed/20.txt']
lab_list = [r'No Restrictions', '2 km', '5 km', '20 km']
lnsty = ['-', '--', '-.', ':']


fig, ax = plt.subplots(1, 1, figsize = (7,4))

for data_path, lab, ls in zip(path_list, lab_list, lnsty):

    data = pd.read_csv(data_path, index_col=False)
    #plt.hist(data.iloc[:, 0], bins = [x for x in range(2, 600)])
    
    c = np.cumsum(np.bincount(data.to_numpy().flatten().astype(np.int)))
    #c = np.bincount(data.to_numpy().flatten().astype(np.int))
    
    data.to_numpy()
    plt.plot(range(1,41), (max(c)-c[:40])/max(c),ls,label = lab)
    #ax.plot(range(len(c)), c,label = lab)
    #plt.hist(data.to_numpy().flatten(), bins = [x for x in range(600)], label = lab)
    
ax.set_xlabel('Distance [km]')
ax.set_ylabel('Complementary cumulative probability function')

#ax.set_xlim(0, 40)
axins = inset_axes(ax, width="40%", height="40%", loc=1)

#plt.ylim(0,200000)
ax.legend(loc=3)

data_path = '../data/processed/1000.txt'
lab = 'No Restrictions'

data = pd.read_csv(data_path, index_col=False)
#plt.hist(data.iloc[:, 0], bins = [x for x in range(2, 600)])

c = np.cumsum(np.bincount(data.to_numpy().flatten().astype(np.int)))

data.to_numpy()
axins.plot((max(c)-c)/max(c),label = lab)
#axins.plot(range(len(c)), c,label = lab)

#plt.hist(data.to_numpy().flatten(), bins = [x for x in range(600)], label = lab)

axins.set_yscale('log')
#axins.set_xscale('log')
axins.set_xlabel('Distance [km]')


#plt.xlabel('Distance [Km]')
#plt.xlim(0, 500)

#plt.legend()
plt.savefig("ccdf.png", dpi = 600)
