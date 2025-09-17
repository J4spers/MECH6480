#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 15:46:22 2025

@author: jasperwetzel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as cm
import time as timer
plt.close('all')


#paramiters
LX = 1.0
LY = 1.0

RHO = 1
MU = 0.01
nu = MU/RHO

#boundaries 

UNORTH = 1.0
USOUTH = 0.0
VEAST = 0.0
VWEST = 0.0

# grid - spacial

NX = 20
NY = 20
dx = LX/NX
dy = LY/NY
dxdy = dx*dy

xnodes = np.linspace(dx/2., LX-dx/2., NX)
ynodes = np.linspace(dy/2., LY-dy/2., NY)
xx, yy = np.meshgrid(xnodes, ynodes, indexing='ij')

#grid - temporal 
DT = 0.001
TOTAL_TIME = 10
num_steps = int(TOTAL_TIME/DT)
PLOT_EVERY = 100

#alocate space 

u = np.zeros((NX+1, NY+2))
ut = np.zeros_like(u)
J_u_x = np.zeros((NX,NY))
J_u_y = np.zeros((NX-1,NY+1))

v = np.zeros((NX+2, NY+1))
vt = np.zeros_like(v)
J_v_x = np.zeros((NX+1,NY-1))
J_v_y = np.zeros((NX,NY))

p = np.zeros((NX+2,NY+2))
prhs = np.zeros_like(p)

uu = np.zeros((NX,NY))
vv = np.zeros((NX,NY))

#initialisation

time=0
u[:,0] = 2.*USOUTH - u[:,1]
u[:,-1] = 2.*UNORTH - u[:,-2]

v[0,:] = 2.*VWEST - v[1,:]
v[-1,:] = 2.*VEAST - v[-2,:]

tic = timer.time()
fig, ax1 = plt.subplots(figsize=(6,6))
for steps in range(num_steps):
    
    #step 1, Calculate u_star and v_star
    
    J_u_x = 0.25*(u[:-1,1:-1] + u[1:,1:-1])**2
    J_u_y = 0.25*(u[1:-1,1:] + u[1:-1,:-1])*\
                 (v[2:-1,:] + v[1:-2,:]) 

    J_u_x -= nu*(u[1:,1:-1]-u[:-1,1:-1])/dx
    J_u_y -= nu*(u[1:-1,1:]-u[1:-1,:-1])/dy
    
    ut[1:-1,1:-1] = u[1:-1,1:-1] - (DT/dxdy)*\
                    (J_u_x[1:,:]*dy - J_u_x[:-1,:]*dy + \
                     J_u_y[:,1:]*dx - J_u_y[:,:-1]*dx)
   
            
    J_v_x = 0.25*(u[:,2:-1]+u[:,1:-2])*(v[1:,1:-1]+v[:-1,1:-1])
    J_v_y = 0.25*(v[1:-1,1:]+v[1:-1,:-1])**2
    
    J_v_x -= nu*(v[1:,1:-1]-v[:-1,1:-1])/dx
    J_v_y -= nu*(v[1:-1,1:]-v[1:-1,:-1])/dy
    
    vt[1:-1,1:-1] = v[1:-1,1:-1] - (DT/dxdy)*\
                    (J_v_x[1:,:]*dy - J_v_x[:-1,:]*dy + \
                     J_v_y[:,1:]*dx - J_v_y[:,:-1]*dx)
            
    ##Boundaries intermediate
    ut[:,0] = 2.*USOUTH - ut[:,1]
    ut[:,-1] = 2.*UNORTH - ut[:,-2]
    vt[0,:] = 2.*VWEST - vt[1,:]
    vt[-1,:] = 2.*VEAST - vt[-2,:]
    

    #step 2, Calulate corrected presure field
    
    divergence = (ut[1:,1:-1]-ut[:-1,1:-1])/dx+\
                 (vt[1:-1,1:]-vt[1:-1,:-1])/dy
                 
    prhs = np.zeros((NX, NY)) 
    prhs[:,:] = RHO * divergence / DT
    for _ in range(50):
        p_next = np.zeros_like(p)
        p_next[1:-1,1:-1] = (-prhs*dxdy**2 + \
                             (p[:-2,1:-1] + p[2:,1:-1])*dy**2 + \
                             (p[1:-1,:-2] + p[1:-1,2:])*dx**2)/\
                             (2.*(dx**2 + dy**2))
        p_next[0,:] = p_next[1,:] 
        p_next[-1,:] = -p_next[-2,:]
        p_next[:,0] = p_next[:,1]
        p_next[:,-1] = p_next[:,-2]
        p = p_next.copy()
    
    #step 3, Correct the velocity field (u and v)
    
    u[1:-1,1:-1] = ut[1:-1,1:-1] - (DT/RHO)*(p[2:-1,1:-1]-\
                                             p[1:-2,1:-1])/dx
    v[1:-1,1:-1] = vt[1:-1,1:-1] - (DT/RHO)*(p[1:-1,2:-1]-\
                                             p[1:-1,1:-2])/dy
    u[:,0] = 2.*USOUTH - u[:,1]
    u[:,-1] = 2.*UNORTH - u[:,-2]
    v[0,:] = 2.*VWEST - v[1,:]
    v[-1,:] = 2.*VEAST - v[-2,:]
    
    time += DT
    if (steps % PLOT_EVERY) == 0:
        uu = 0.5*(u[:-1,1:-1]+u[1:,1:-1])
        vv = 0.5*(v[1:-1,:-1]+v[1:-1,1:])
        ax1.clear()
        #vel_fig = ax1.pcolormesh(xx,yy,np.sqrt(uu**2+vv**2), shading='auto', cmap='jet')
        vel_fig = ax1.contourf(xx,yy, np.sqrt(uu**2+vv**2), levels=20, cmap='jet')
        vel_bar = plt.colorbar(vel_fig, ax=ax1)
        ax1.quiver(xx,yy,uu,vv)
        ax1.set_xlim(0,LX)
        ax1.set_ylim(0,LY)
        ax1.set_aspect('equal')
        plt.pause(0.01)
        vel_bar.remove()
        
plt.show()