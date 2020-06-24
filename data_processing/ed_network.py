#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 06:37:21 2020

@author: roryh
"""

import csv
import pandas as pd
import datetime
import numpy as np
import matplotlib.pyplot as plt

# hold info on every Electoral Division (ed)
ed_basic_info = pd.read_csv('../data/raw/ED_Basic_Info.csv')

ed_names = ed_basic_info['Electoral Division'].to_numpy()
# maps each ed to it's index in ed_basic_info
ed_id_map = {ed:x for ed, x in zip(ed_names, ed_basic_info.index)}

# has the number of commuter between each ed
ed_links = pd.read_csv('../data/raw/ED_Used_Link_Info.csv')

# a dense matrix to hold the number of cummters from ed i to j
travel_mat = np.zeros((ed_names.size, ed_names.size))
for row in ed_links.to_numpy():
    src = ed_id_map[row[0]]
    
    if row[2] == 'No fixed place of work':
        continue
    if row[2] == 'Work/school from home':
        continue
    
    dst = ed_id_map[row[2]]
    w = int(row[4])
    travel_mat[src, dst] += w

# populate a data frame that holds the proprotion of commuters in i that travel to j
# so the weight is just the probability of someone in i travelling to j
travel_dict = {'src':[], 'dst':[], 'weight':[]}
for i, row in enumerate(travel_mat):
    w_sum = np.sum(row)
    for j, val in enumerate(row):
        if val == 0:
            continue
        travel_dict['src'].append(i)
        travel_dict['dst'].append(j)
        travel_dict['weight'].append(val)
        
vertices_df = pd.DataFrame(ed_basic_info)
vertices_df = vertices_df.drop(['Electoral Division', 'COUNTY', 'ED Area'], axis=1)
vertices_df . columns = ['long', 'lat', 'population', 'commuters']
vertices_df.to_csv('../data/processed/ed_vertices.csv', index = True)

edges_df = pd.DataFrame(travel_dict)
edges_df.to_csv('../data/processed/ed_edges.csv', index = False)
