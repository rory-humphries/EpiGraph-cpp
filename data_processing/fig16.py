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
    c = np.cumsum(np.bincount(data.to_numpy().flatten().astype(np.int)))
    data.to_numpy()
    plt.plot(range(1,41), (max(c)-c[:40])/max(c),ls,label = lab)
    
ax.set_xlabel('Distance [km]')
ax.set_ylabel('Complementary cumulative probability function')
axins = inset_axes(ax, width="40%", height="40%", loc=1)
ax.legend(loc=3)
data_path = '../data/processed/1000.txt'
lab = 'No Restrictions'
data = pd.read_csv(data_path, index_col=False)
c = np.cumsum(np.bincount(data.to_numpy().flatten().astype(np.int)))
data.to_numpy()
axins.plot((max(c)-c)/max(c),label = lab)
axins.set_yscale('log')
axins.set_xlabel('Distance [km]')
plt.savefig("output", dpi = 600)
