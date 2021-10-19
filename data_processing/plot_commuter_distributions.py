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

data_path = '../data/processed/1000.txt'

output_name = 'distance_distribution'
fig_name = output_name + '.png'
data_name = output_name + '.csv'
#####################

path_list = ['../data/processed/trav_mat_1000.csv',
             '../data/processed/trav_mat_2.csv',
             '../data/processed/trav_mat_5.csv',
             '../data/processed/trav_mat_20.csv']
lab_list = [r'No restrictions', '2 km', '5 km', '20 km']
lnsty = ['-', '--', '-.', ':']

fig, ax = plt.subplots(1, 1, figsize = (7,4))

for data_path, lab, ls in zip(path_list, lab_list, lnsty):
    

    #data = pd.read_csv(data_path, index_col=False)
    data = np.genfromtxt(data_path, delimiter=',')[:,:-1]
    #plt.hist(data.iloc[:, 0], bins = [x for x in range(2, 600)])
    
    #c = np.cumsum(np.bincount(data.flatten().astype(np.int)))
    c = np.bincount(data.flatten().astype(np.int))
    
    #plt.plot((c)/max(c),label = lab)
    ax.plot(range(1, len(c)+1), c,ls, label = lab)
    #plt.hist(data.flatten(), bins = [x for x in range(600)], label = lab)
#plt.xlim(0,40)
#plt.yscale('log')
    #plt.show()
    
ax.set_xlabel(r'$N_{ij}$')
#ax.set_xlim(0, 40)
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xscale('log')

# axins = inset_axes(ax, width="40%", height="40%", loc=1)

# #plt.ylim(0,200000)
# ax.legend(loc=2)

# for data_path, lab in zip(path_list, lab_list):

#     data = pd.read_csv(data_path, index_col=False)
#     #plt.hist(data.iloc[:, 0], bins = [x for x in range(2, 600)])
    
#     c = np.cumsum(np.bincount(data.to_numpy().flatten().astype(np.int)))
    
#     data.to_numpy()
#     axins.plot((c)/max(c),label = lab)
#     #axins.plot(range(len(c)), c,label = lab)

#     #plt.hist(data.to_numpy().flatten(), bins = [x for x in range(600)], label = lab)

# axins.set_yscale('log')
# axins.set_xscale('log')

# #plt.xlabel('Distance [Km]')
# #plt.xlim(0, 500)

plt.legend()
plt.savefig("cumulative_commuter_dists.png", dpi = 600)
