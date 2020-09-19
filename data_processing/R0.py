#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 11:11:47 2020

@author: roryh
"""
import os

direc = '../mats/'
files = []
for filename in os.listdir(direc):
    files += [filename]
    
edges = pd.read_csv('../mats/0.csv')
vertices = pd.read_csv('../output/0.csv')

S = vertices.S.to_numpy()
I = vertices.I.to_numpy()
N = vertices.N.to_numpy()


import scipy.sparse as spr

src = df.source.to_numpy()
dst = df.destination.to_numpy()
pop = df.population.to_numpy()

mat = spr.coo_matrix((pop, (src, dst)))

mat2 = mat.copy()

mat2 = mat2.tocsr()

for i in range(mat.shape[0]):
    nout = sum(mat.getrow(i).data)
    mat2[i].data -= nout
    
mat2 = mat2.tocsc()

for i in range(mat.shape[1]):
    nin = sum(mat.getcol(i).data)
    mat2[:,i].data -= nin

sum(mat[0].data)


for row, col, dat in (mat.row, mat.col, mat.data):
    N[row] -= dat
    N[col] += dat
