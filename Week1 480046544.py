#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 09:53:38 2025

@author: jasperwetzel
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import git
from datetime import datetime
from git import Repo

## Importing data
data = pd.read_csv("curve.data", sep='  ', names=['x', 'y'], header=0)
x = data['x']
y = data['y']


## Plotting datta and timestamp

Time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

data.plot(kind = 'line', x = 'x', y = 'y', grid = 'True', title = 'Example plot of data in curve.data')
plt.annotate(Time, xy=(0.7,0.95), xycoords='figure fraction', annotation_clip=False)


## Git repo
repo = Repo('.', search_parent_directories=True)
revsha = repo.head.object.hexsha[:8]
plt.annotate(f"[rev {revsha}]", xy=(0.05,0.95), xycoords='figure fraction', annotation_clip=False)

plt.show()

Fwrd_dy = []
def Forward_diff(data):
    for n in range(len(data)-1):
        fl = (y[n+1] - y[n]) / (x[n+1] - x[n])
        Fwrd_dy.append(fl)

Central_dy = []
def Central_diff(data):
    for n in range(1,len(data)-1):
        fl = (y[n+1] - y[n-1]) / 2/(x[n+1] - x[n])
        Central_dy.append(fl)
        
Forward_diff(y)
Central_diff(y)

plt.plot(x[: len(y) - 1], Fwrd_dy, 
         color='red', label='Forward Difference')

plt.plot(x[: len(y) - 2], Central_dy, 
         color='green', label='Central Difference')
plt.plot(x, -1* np.sin(x),color='blue', label='Actual' )
plt.grid()
plt.legend()

repo = Repo('.', search_parent_directories=True)
revsha = repo.head.object.hexsha[:8]
plt.annotate(f"[rev {revsha}]", xy=(0.05,0.95), xycoords='figure fraction', annotation_clip=False)

plt.annotate(Time, xy=(0.7,0.95), xycoords='figure fraction', annotation_clip=False)

plt.show()

a = -1* np.sin(x[:99]) - Fwrd_dy
b = -1* np.sin(x[:98]) - Central_dy

plt.plot(x[: len(a)], a, 
         color='red', label='Forward error')

plt.plot(x[: len(b)], b, 
         color='green', label='Central error')

plt.grid()
plt.legend()

repo = Repo('.', search_parent_directories=True)
revsha = repo.head.object.hexsha[:8]
plt.annotate(f"[rev {revsha}]", xy=(0.05,0.95), xycoords='figure fraction', annotation_clip=False)

plt.annotate(Time, xy=(0.7,0.95), xycoords='figure fraction', annotation_clip=False)

plt.show()