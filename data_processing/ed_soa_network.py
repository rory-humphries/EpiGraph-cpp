#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 21:05:36 2020

@author: roryh
"""

import csv
import pandas as pd
import datetime
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd

from pyproj import Proj, transform

"""
This script outputs all the data needed to run network models on the combined
ROI electoral division data and the NI super output area data.
"""


"""
Output the vertex id's and there location in long/lat
"""

# hold info on every Electoral Division (ed) and super output area (soa)
ed_soa_df = pd.read_csv('../data/raw/Joined_Pop_Data_CSO_NISRA.csv')

# maps each ed to and soa to it's index in ed_soa_df
ed_soa_id_map = {ed:x for ed, x in zip(ed_soa_df['Electoral Division'], ed_soa_df.index)}

soa_counties = pd.read_csv('../data/raw/Geographic Data (statistical geographies).csv')
soa_counties.index = soa_counties['SOA Code']

ed_counties = pd.read_csv('../data/raw/ED_Basic_Info.csv')
ed_counties.index = ed_counties['Electoral Division'].str.split(expand=True)[0]

county_set = []
county_set.extend(np.unique(ed_counties['COUNTY'].to_numpy()))
county_dict = {x:y for x,y in zip(county_set, range(len(county_set)))}

ed_soa_df.index = ed_soa_df['Electoral Division'].str.split(expand=True)[0]

ed_soa_df['type'] = len(county_dict)
for e in ed_counties.index:
    ed_soa_df['type'][e] = county_dict[ed_counties['COUNTY'][e]]
    

vertices_df = pd.DataFrame(ed_soa_df)
vertices_df = vertices_df.drop(['Electoral Division'], axis=1)
vertices_df . columns = ['long', 'lat', 'population', 'type']

from pyproj import Transformer
transformer = Transformer.from_crs(2157, 4326, always_xy = True)

long_list = []
lat_list = []
for pt in transformer.itransform(zip(vertices_df.long, vertices_df.lat)):
    long_list.append(pt[0])
    lat_list.append(pt[1])

vertices_df.long = long_list
vertices_df.lat = lat_list

vertices_df.to_csv('../data/processed/ed_soa_vertices.csv', index=True)

vertices_df[["long", "lat"]].to_csv('../data/processed/ed_soa_long_lat.csv', index=False)
vertices_df[["population"]].to_csv('../data/processed/ed_soa_population.csv', index=False)
vertices_df[["type"]].to_csv('../data/processed/ed_soa_county.csv', index=False)

"""
Output the probability distribution of travel distances from the ed edges
"""

# read in electoral division data and commuting data
ed_df = pd.read_csv('../data/raw/ED_Basic_Info.csv')
ed_names = ed_df['Electoral Division'].to_numpy()

# maps each ed to it's index in ed_basic_info
ed_id_map = {ed:x for ed, x in zip(ed_df['Electoral Division'], ed_df.index)}

# has the number of commuter between each ed
ed_links = pd.read_csv('../data/raw/ED_Used_Link_Info.csv')

# will hold the probability of a commuter travelling a certain (rounded to int) distance
ed_dist = dict()
for i,j,k in zip(ed_links['Electoral Division'].to_numpy(), ed_links.Distance.to_numpy(), ed_links['No. of Commuters'].to_numpy()):
    if int(j) not in ed_dist.keys():
        ed_dist[int(j)] = int(k)
    else:
        ed_dist[int(j)] += int(k)  
tot = sum(ed_dist.values())
for j in ed_dist.keys():
    ed_dist[j]/=tot

# insert the data into a dataframe and write to file
dist = list()
probs = list()

for i,j in ed_dist.items():
    dist.append(i)
    probs.append(j)
    
distance_distribution = pd.DataFrame({'value':dist, 'probability':probs})

distance_distribution.to_csv('../data/processed/ed_distance_probs.csv', index = False)

"""
Output the probability distribution of commuter proportions
"""
bins = np.linspace(0,1,50)
bin_indices = np.digitize(ed_df['No. of Commuters'].to_numpy()/ed_df.Population.to_numpy(), bins)
bin_count = np.bincount(bin_indices)
probs = bin_count/np.sum(bin_count)

commuter_distribution = pd.DataFrame({'value':bins[:len(probs)], 'probability':probs})
commuter_distribution.to_csv('../data/processed/ed_commuter_probs.csv', index = False)


"""
Output the matrix which gives the probaility of a vertex travelling to another based 
on the dostribution of travel distances from electoral divisions
"""

bands = pd.read_csv("../data/raw/Distance_Distribution_by_Ratio.csv")

bands_dic = {x.split('-')[-1].strip():np.zeros(550) for x in bands["Ratio Band"]}
bands_dic = {x.split('-')[-1].strip():np.zeros(550) for x in bands["Ratio Band"]}

for b, d, p in zip(bands["Ratio Band"], bands["Distance"], bands["Prob"]):
    bands_dic[b.split('-')[-1].strip()][int(d)] = p
    
for k in bands_dic.keys():
    prev = 0
    for i, val in enumerate(bands_dic[k]):
        if val == 0:
            bands_dic[k][i] = prev
        else:
            prev = val

ratios = np.array([0.1, 1, 3800])
ratios_str = np.array(['0.1', '1', '3800'])
pop_list = vertices_df.population.to_numpy()
vertex_travel_mat = np.zeros((len(vertices_df.index), len(vertices_df.index)))

for src in range(len(vertices_df.index)):
    print(src)
    for dst in range(len(vertices_df.index)):
        if src == dst:
            continue
        else:
            ratio = pop_list[src]/pop_list[dst]
            ratio = ratios_str[np.searchsorted(ratios, ratio, 'left')]
    
            d = int(distance_haversine(lat_list[src], long_list[src], 
                                       lat_list[dst], long_list[dst])/1000.0)
            
            #d = bands_dic[ratio][d]
            
            vertex_travel_mat[src][dst] = bands_dic[ratio][d]

np.savetxt('../data/processed/ed_soa_travel_prob_mat_ratio_bands.csv', vertex_travel_mat, delimiter=',', fmt = '%f')

def radians(x):
    return x*np.pi/180.0

def distance_haversine(lat1, lon1, lat2, lon2):
    # approximate radius of earth in km
    R = 6373.0 * 1000

    lat1 = radians(lat1)
    lon1 = radians(lon1)
    lat2 = radians(lat2)
    lon2 = radians(lon2)

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    distance = R * c

    return distance

dist_probs_map = {x:y for x,y in zip(distance_distribution.distance, distance_distribution.probability)}
for u in range(500):
    if u not in dist_probs_map.keys():
        dist_probs_map[u]=0
        
vertex_travel_mat = np.zeros((len(vertices_df.index), len(vertices_df.index)))

for src in vertices_df.index:
    #print(src)
    tot = 0
    for dst in vertices_df.index:
        if src == dst:
            continue
        else:
            tot += dist_probs_map[int(distance_haversine(lat_list[src], 
                                                             long_list[src], 
                                                             lat_list[dst], 
                                                             long_list[dst])/1000.0)]
            vertex_travel_mat[src][dst] = tot
            
for x in vertices_df.index:
    vertex_travel_mat[x] /= np.sum(vertex_travel_mat[x])
                                   
np.savetxt('../data/processed/ed_soa_travel_prob_mat.csv', vertex_travel_mat, delimiter=',', fmt = '%f')
            
