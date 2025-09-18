#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 13 20:45:45 2025

@author: jasperwetzel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as cm
import time as timer
from datetime import datetime
from git import Repo
import time as timer
plt.close('all')

Time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

#plt.annotate(Time, xy=(0.7,0.95), xycoords='figure fraction', annotation_clip=False)

## Git repo

repo = Repo('.', search_parent_directories=True)
revsha = repo.head.object.hexsha[:8]
#plt.annotate(f"[rev {revsha}]", xy=(0.05,0.95), xycoords='figure fraction', annotation_clip=False)

#paramiters
LX = 0.1
LY = 0.1

RHO = 1.0
NU = 1.0e-5
mu = NU*RHO
PR = 0.71
RA = [1e5, 1e6]
TNORTH = 300
TSOUTH = 400
TIN = 350
UIN = 0
VIN = 0
G = [0,10]
H = 0.1
alpha = NU/PR
beta_max = alpha * RA[1] * NU / (H**3 * -G[1] * (TNORTH - TSOUTH))
beta_min = alpha * RA[0] * NU / (H**3 * -G[1] * (TNORTH - TSOUTH))

#boundaries 

def Nu(flux):
    return H/(TSOUTH-TNORTH) * flux[:,0]




def noslip_u(vel):
    vel[:,-1] = -vel[:,-2]
    vel[:,0] = -vel[:,1]
    vel[0,:] = 0.0
    vel[-1,:] = 0.0
    return

def noslip_v(vel):
    vel[:,0] = 0.0
    vel[:,-1] = 0.0
    vel[-1,:] = -vel[-2,:]
    vel[0,:] = -vel[1,:]
    return

def T_north(temp):
    temp[:,-1] = TNORTH*2 - temp[:,-2]
    return

def T_south(temp):
    temp[:,0] = TSOUTH*2 - temp[:,1]
    return

def T_east(temp):
    temp[-1,1:-1] = -temp[-2,1:-1]
    return

def T_west(temp):
    temp[0,1:-1] = -temp[-1,1:-1]
    return

# grid - spacial

NX = 40
NY = 40
dx = LX/NX
dy = LY/NY
dxdy = dx*dy

xnodes = np.linspace(dx/2., LX-dx/2., NX)
ynodes = np.linspace(dy/2., LY-dy/2., NY)
xx, yy = np.meshgrid(xnodes, ynodes, indexing='ij')

#grid - temporal
DT = 1e-5
TOTAL_TIME = 10#*DT
num_steps = int(TOTAL_TIME/DT)
PLOT_EVERY = 100

#alocate space 

T = np.ones((NX+2,NY+2))*TIN

u = np.zeros((NX+1, NY+2))
ut = np.zeros_like(u)
J_u_x = np.zeros((NX,NY))
J_u_y = np.zeros((NX-1,NY+1))

v = np.zeros((NX+2, NY+1))
vt = np.zeros_like(v)
J_v_x = np.zeros((NX+1,NY-1))
J_v_y = np.zeros((NX,NY))

T_NS = np.zeros((NX,NY+1))
T_NS_c = np.zeros_like(T_NS)
T_EW = np.zeros((NX+1,NY))
T_EW_c = np.zeros_like(T_EW)

p = np.zeros((NX+2,NY+2))
prhs = np.zeros_like(p)

uu = np.zeros((NX,NY))
vv = np.zeros((NX,NY))

#initialisation

time=0

noslip_u(u)
noslip_v(v)
T_north(T)
T_south(T)
T_east(T)
T_west(T)

tic = timer.time()
fig, ax1 = plt.subplots(figsize=(6,6))
for steps in range(num_steps):
    
    #step 0, Calculate tempratures
    
    #T[1:-1,1:-1] = T[1:-1,1:-1] - (DT/dx) * (u[1:,1:-1]*(T[1:-1,1:-1]+T[2:,1:-1]) - u[:-1,1:-1]*(T[1:-1,1:-1]+T[:-2,1:-1])) \
                                #- ((DT/dy) * (v[1:-1,1:]*(T[1:-1,1:-1]+T[1:-1,2:]) - v[1:-1,:-1]*(T[1:-1,1:-1]+T[1:-1,:-2]))) \
                                    #+ (alpha * DT * ((T[2:,1:-1]   - 2*T[1:-1,1:-1] + T[:-2,1:-1])   / dx**2
                                                   #+ (T[1:-1,2:]  - 2*T[1:-1,1:-1] + T[1:-1,:-2])   / dy**2)) 
    
    #step 1, Calculate u_star and v_star
    
    J_u_x = 0.25*(u[:-1,1:-1] + u[1:,1:-1])**2
    J_u_y = 0.25*(u[1:-1,1:] + u[1:-1,:-1])*\
                 (v[2:-1,:] + v[1:-2,:]) 

    J_u_x -= NU*(u[1:,1:-1]-u[:-1,1:-1])/dx
    J_u_y -= NU*(u[1:-1,1:]-u[1:-1,:-1])/dy
    
    ut[1:-1,1:-1] = u[1:-1,1:-1] - (DT/dxdy) *\
                    (J_u_x[1:,:]*dy - J_u_x[:-1,:]*dy + \
                     J_u_y[:,1:]*dx - J_u_y[:,:-1]*dx) #+ 1*DT
   
            
    J_v_x = 0.25*(u[:,2:-1]+u[:,1:-2])*(v[1:,1:-1]+v[:-1,1:-1])
    J_v_y = 0.25*(v[1:-1,1:]+v[1:-1,:-1])**2
    
    J_v_x -= NU*(v[1:,1:-1]-v[:-1,1:-1])/dx
    J_v_y -= NU*(v[1:-1,1:]-v[1:-1,:-1])/dy
    
    vt[1:-1,1:-1] = v[1:-1,1:-1] - (DT/dxdy)*\
                    (J_v_x[1:,:]*dy - J_v_x[:-1,:]*dy + \
                     J_v_y[:,1:]*dx - J_v_y[:,:-1]*dx) + (1 - (G[1] * beta_max*(TIN-(0.5*(T[1:-1,1:-2]+T[1:-1,2:-1])))))*DT
    
     
    #T[1:-1,1:-1] = T[1:-1,1:-1] + alpha * DT * ((T[2:,1:-1]   - 2*T[1:-1,1:-1] + T[:-2,1:-1])   / dx**2
                                                    #+ (T[1:-1,2:]  - 2*T[1:-1,1:-1] + T[1:-1,:-2])   / dy**2)
    ##Boundaries intermediate
    noslip_u(ut)
    noslip_v(vt)
    

    #step 2, Calulate corrected presure field
    
    divergence = (ut[1:,1:-1]-ut[:-1,1:-1])/dx+\
                 (vt[1:-1,1:]-vt[1:-1,:-1])/dy
                 
    prhs = np.zeros((NX, NY)) 
    prhs[:,:] = RHO * divergence / DT
    prhs -= prhs.mean()
    for _ in range(200):
        p_next = np.zeros_like(p)
        p_next[1:-1,1:-1] = (-prhs*dxdy**2 + \
                             (p[:-2,1:-1] + p[2:,1:-1])*dy**2 + \
                             (p[1:-1,:-2] + p[1:-1,2:])*dx**2)/\
                             (2.*(dx**2 + dy**2))
        p_next[0,:] = -p_next[1,:] 
        p_next[-1,:] = -p_next[-2,:]
        p_next[:,0] = -p_next[:,1]
        p_next[:,-1] = -p_next[:,-2]
        #p_next[1,1] = 0.0
        p = p_next.copy()
    
    #step 3, Correct the velocity field (u and v)
    
    u[1:-1,1:-1] = ut[1:-1,1:-1] - (DT/RHO)*(p[2:-1,1:-1]-\
                                             p[1:-2,1:-1])/dx
    v[1:-1,1:-1] = vt[1:-1,1:-1] - (DT/RHO)*(p[1:-1,2:-1]-\
                                             p[1:-1,1:-2])/dy
    noslip_u(u)
    noslip_v(v)
    
    # Diffusion
    T_EW[1:-1,:] = alpha * (T[2:-1,1:-1] - T[1:-2,1:-1]) /dx
    T_NS[:,1:-1] = alpha * (T[1:-1,2:-1] - T[1:-1, 1:-2]) / dy
    
    # Convection
    T_EW_c[1:-1,:] = (T[2:-1,1:-1] + T[1:-2,1:-1]) / 2
    T_NS_c[:,1:-1] = (T[1:-1,2:-1] + T[1:-1,1:-2]) / 2
    
    # Boundary 
    # top wall
    T_NS[:,-1] = alpha * (TNORTH - T[1:-1, -2]) / (dy/2)
    T_NS_c[:,-1] = TNORTH
    # bottom wall
    T_NS[:,0] = alpha * (T[1:-1, 1] - TSOUTH) / (dy/2)
    T_NS_c[:,0] = TSOUTH
    # east wall
    T_EW[-1,:] = 0
    T_e = T[-3,1:-1] + 1.5 * dx * ((T[-3,1:-1]-T[-2,1:-1])/dx)
    T_EW_c[-1,:] = T_e
    # west wall
    T_EW[0,:] = 0
    T_w = T[2,1:-1] + 1.5 * dx * ((T[1,1:-1]-T[2,1:-1])/dx)
    T_EW_c[0,:] = T_w
    
    T[1:-1,1:-1] = T[1:-1,1:-1] + DT * (1/(dx*dy)) * \
        (dy * T_EW[1:,:] - dy * T_EW[:-1,:] + dx * T_NS[:,1:] - dx * T_NS[:,:-1]) -\
            DT * (1/(dx*dy)) * (T_EW_c[1:,:] * u[1:,1:-1] * dy - T_EW_c[:-1,:] * u[:-1,1:-1] * dy +\
                                T_NS_c[:,1:] * v[1:-1,1:] * dx - T_NS_c[:,:-1] * v[1:-1,:-1] * dx)
                
    #T_east(T)
    #T_north(T)
    #T_south(T)
    #T_west(T)
    '''n 
    T_EW[:,:] = alpha * (T[1:,1:-1] - T[:-1,1:-1])/dx
    T_NS[:,:] = alpha * (T[1:-1,1:] - T[1:-1,:-1])/dy
    
    T_EW_c[1:-1,:] = 0.5*((T[2:-1,1:-1] + T[1:-2,1:-1]) * u[1:-1,1:-1] -\
                  (T[1:-2,1:-1] + T[2:-1,1:-1]) * u[1:-1,1:-1])*dy
        
    T_NS_c[:,1:-1] = 0.5*((T[1:-1,2:-1] + T[1:-1,1:-2]) * v[1:-1,1:-1] -\
                  (T[1:-1,1:-2] + T[1:-1,2:-1]) * v[1:-1,1:-1])*dx
    
    #Boundaries    
    T_NS_c[:,0] = TSOUTH* v[1:-1,0] /dy
    T_NS_c[:,-1] = TNORTH* v[1:-1,-1] /dy
    
    T_EW_c[0,:] = u[0,1:-1] * (T[0,1:-1] + T[1,1:-1])*0.5 /dx
    T_EW_c[-1,:] = u[-1,1:-1] * (T[-1,1:-1] + T[-2,1:-1])*0.5 /dx
    
    T[1:-1,1:-1] = T[1:-1,1:-1] - DT*(((T_EW_c[1:,:]-T_EW_c[:-1,:]) +\
                                       (T_NS_c[:,1:]-T_NS_c[:,:-1])) +\
                                      ((T_EW[1:,:]-T_EW[:-1,:]) +\
                                       (T_NS[:,1:]-T_NS[:,:-1]))) /dxdy
    
    T_north(T)
    T_south(T)
    T_east(T)
    T_west(T)
   '''
    
    time += DT
    if (steps % PLOT_EVERY) == 0:
        
        uu = 0.5*(u[:-1,1:-1]+u[1:,1:-1])
        vv = 0.5*(v[1:-1,:-1]+v[1:-1,1:])
        ax1.clear()
        #vel_fig = ax1.pcolormesh(xx,yy,np.sqrt(uu**2+vv**2), shading='auto', cmap='jet')
        T_rounding = np.around(T[1:-1,1:-1],2)
        vel_fig = ax1.contourf(xx,yy, T[1:-1,1:-1], levels=20, cmap='hot',label=f't={time:.2f}s')
        vel_bar = plt.colorbar(vel_fig, ax=ax1)
        plt.annotate(Time, xy=(0.7,0.95), xycoords='figure fraction', annotation_clip=False)
        plt.annotate(f't={time:.3f}s', xy=(0.4,0.85), xycoords='figure fraction', annotation_clip=False)

        ax1.quiver(xx,yy,uu,vv)
        ax1.set_xlim(0,LX)
        ax1.set_ylim(0,LY)
        ax1.set_aspect('equal')
        plt.pause(0.01)
        vel_bar.remove()
        
plt.show()