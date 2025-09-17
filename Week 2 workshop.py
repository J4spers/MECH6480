#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 17:10:15 2025

@author: jasperwetzel
"""

import numpy as np
import matplotlib.pyplot as plt

# Heat diffusion through a bar

def Steady_Diffusion(k, dx, num_cells, T_A, T_B):
    a_W = k /dx
    a_E = k/dx
    a_P = a_E + a_W
    
    S_u = -2.*k / dx
    a_p_star = 3*k /dx
    
    A = np.diag(np.repeat(-a_P, num_cells)) + \
        np.diag(np.repeat(a_E, num_cells-1), 1) + \
        np.diag(np.repeat(a_W, num_cells-1), -1)
    A[0,0] = -1*a_p_star
    A[-1,-1] = -1*a_p_star
    b = np.zeros(num_cells)
    b[0] = T_A * S_u
    b[-1] = T_B * S_u
    return A, b.transpose()

def Steady_Convection(F, dx, num_cells, phi_A, phi_B):
    a_W = -F/2.
    a_E = F/2.
    a_P = a_E + a_W
    S_u = F
    a_p_star = F/2.
    A = np.diag(np.repeat(-a_P, num_cells)) + \
        np.diag(np.repeat(a_E, num_cells-1), 1) + \
        np.diag(np.repeat(a_W, num_cells-1), -1)
    A[0,0] = 1*a_p_star
    A[-1,-1] = -1*a_p_star
    b = np.zeros(num_cells)
    b[0] = S_u * phi_A
    b[-1] = -S_u * phi_B
    return A, b.transpose()

# Values 
k = 1000
length = 0.5
T_A = 100
T_B = 500

# Grid Gen
num_cells = 5
dx = length / num_cells
x_locations = np.linspace(dx/2., length-dx/2.,num_cells)
A, b = Steady_Diffusion(k, dx, num_cells, T_A, T_B)
T = np.linalg.solve(A, b)

#Plotting 
plt.plot(x_locations, T, 'bo')
plt.plot([0, length], [T_A, T_B], 'r+', label = 'boundary')
plt.grid()
plt.show()

# Values 2
gamma =0.1
length = 1
phi_A = 1
phi_B = 0
roh = 1
u = 0.1

num_cells = 5
dx = length / num_cells
x_locations = np.linspace(dx/2., length-dx/2.,num_cells)
D, b_D = Steady_Diffusion(gamma, dx, num_cells, phi_A, phi_B)
A, b_A = Steady_Convection(roh*u, dx, num_cells, phi_A, phi_B)
T = np.linalg.solve(D-A, b_D - b_A)

plt.plot(x_locations, T, 'bo')
plt.plot([0, length], [phi_A, phi_B], 'r+', label = 'boundary')
plt.grid()
plt.show()