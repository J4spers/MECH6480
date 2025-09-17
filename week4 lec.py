#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 15:48:45 2025

@author: jasperwetzel
"""

import numpy as np
import matplotlib.pyplot as plt

#paarmaters 
L = 1.5
u = 2.0
p= 1.0
G = 0.03

def s(x):
    if 0 <= x <= 0.6:
        return -200 * x +100
    elif 0.6 < x <= 0.8:
        return 100 * x -80
    else:
        return 0

#grid space
num_cells = 20
dx = L/ num_cells
x_locs = np.linspace(0.5*dx, L-0.5*dx, num_cells)

#gris time
dt = 0.01
total_time = 2
num_steps = int(total_time/dt)

#initial conditions 
concentration = np.zeros(num_cells)
time = 0 

#fluxes
F = p*u
D = G /dx

# matrix 
B = np.diag(np.repeat(p*dx/dt, num_cells))
A = np.diag(np.repeat(p*dx/dt - F- 2*D, num_cells)) + \
    np.diag(np.repeat(D, num_cells-1),1) + \
    np.diag(np.repeat(D+F, num_cells-1),-1)
A[0,0] = p*dx/dt - F -3*D
A[-1,-1] = p*dx/dt - F -D

b = np.array([s(x) * dx for x in x_locs])

update_matrix = np.linalg.inv(B) @ A # or use np.matmul()
source_vector = np.linalg.inv(B) @ b.T

fig, ax = plt.subplots(figsize=(10,6))
#fig, ax = plt.subplot()
for n in range(num_steps):
    concentration = update_matrix @ concentration + source_vector
    time += dt
    ax.cla()
    ax.plot(x_locs, concentration, label=f't={time:.2f}s')
    ax.plot([0,L],[0,concentration[-1]], 'ko', label = 'boundary con')
    ax.grid()
    ax.legend()
    #plt.pause(0.0001)
plt.show()


    