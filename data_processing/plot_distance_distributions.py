#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 11:21:40 2020

@author: roryh
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
######################
save_figs = True

data_path = '../data/processed/1.txt'

output_name = 'distance_distribution'
fig_name = output_name + '.png'
data_name = output_name + '.csv'
#####################

data = pd.read_csv(data_path, index_col=False)
plt.hist(data.iloc[:, 0], bins = [x for x in range(600)])

plt.xlim(0, 10)
#plt.yscale('log')
