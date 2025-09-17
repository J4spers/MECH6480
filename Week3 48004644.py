#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 09:31:22 2025

@author: jasperwetzel
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
from git import Repo

#Making matrix A and solution vector b

def Poisson(k, dr, num_cells, dp):
    r_P = -2
    S_u = (1/k) * dp * dr * dr

    main_diag  = np.arange(1, num_cells+1) * r_P
    upper_diag = np.arange(1, num_cells)
    lower_diag = np.arange(1, num_cells)

    A = (
        np.diag(main_diag + 1) +
        np.diag(upper_diag, 1) +
        np.diag(lower_diag, -1)
    )

# Boundary Conditions
    A[0,0] = -1

    b = S_u * (np.arange(num_cells)+ dr/2).reshape(-1, 1)
    b[0] = S_u
    return A, b

#Values for p

radius = 0.1
num_cells = 20
dr = radius / num_cells
k = 1e-3
dp = -0.02

r_locations = np.linspace(dr/2., radius-dr/2.,num_cells)
A, b = Poisson(k, dr, num_cells, dp)
T = np.linalg.solve(A, b)

#Anylitical solution

U = (dp/(4*k) * (r_locations**2 - radius**2))

# Plotting and Time-Stamp

Time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

plt.plot(r_locations, T, 'bo', label = 'Aproximate Solution')
plt.xlabel('Radius [m]')
plt.ylabel('Velocity [m/s]')
plt.plot([0, radius], [T[0], T[-1]], 'r+', label = 'Boundary')
plt.grid()

plt.annotate(Time, xy=(0.7,0.95), xycoords='figure fraction', annotation_clip=False)
repo = Repo('.', search_parent_directories=True)
revsha = repo.head.object.hexsha[:8]
plt.annotate(f"[rev {revsha}]", xy=(0.05,0.95), xycoords='figure fraction', annotation_clip=False)

plt.plot(r_locations, U, label = 'Analytical Soution', color = 'r')
plt.legend()

plt.show()

