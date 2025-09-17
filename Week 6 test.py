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
import psutil, os
plt.close('all')

udata = np.genfromtxt('exper-xvel-re230-xh0.xy', names=True, delimiter=' ')
vdata = np.genfromtxt('exper-xvel-re230-xh6.xy', names=True, delimiter=' ')


#paramiters
LX = 1.0
LY = 0.25

RHO = 1
MU = 0.01
G = 0.1
nu = MU/RHO

#boundaries 

UNORTH = 1.0
USOUTH = 0.0
VEAST = 0.0
VWEST = 0.0

# grid - spacial

NX = 80
NY = 20
dx = LX/NX
dy = LY/NY
dxdy = dx*dy

xnodes = np.linspace(dx/2., LX-dx/2., NX)
ynodes = np.linspace(dy/2., LY-dy/2., NY)
xx, yy = np.meshgrid(xnodes, ynodes, indexing='ij')

#grid - temporal 
DT = 0.00001
#TOTAL_TIME = 100
#num_steps = int(TOTAL_TIME/DT)
NUM_STEPS = 5000
PLOT_EVERY = 100
N_PRESSURE_ITS = 400
UIN = 0

#alocate space 

u = np.ones((NX+1, NY+2))*UIN
ut = np.zeros_like(u)
J_u_x = np.zeros((NX,NY))
J_u_y = np.zeros((NX,NY+1))

v = np.zeros((NX+2, NY+1))
vt = np.zeros_like(v)
J_v_x = np.zeros((NX+1,NY-1))
J_v_y = np.zeros((NX,NY))

p = np.zeros((NX+2,NY+2))
prhs = np.zeros_like(p)
p_next = np.zeros_like(p)

uu = np.zeros((NX,NY))
vv = np.zeros((NX,NY))

#Boundary conditions

def inflow_outflow_xmomentum_noscale(vel):
    #vel[0,1:-1] = UIN
   # u_in =np.sum(vel[0,1:-1])
    #u_out = np.sum(vel[-2,1:-1])
    #vel[-1,1:-1] = vel[-2,1:-1]#*u_in/u_out
    vel[0,:] = vel[-1,:]
    return

def inflow_outflow_xmomentum_noscale(vel):
    vel[0,1:-1] =  UIN #vel[-1,1:-1]
    u_in =np.sum(vel[0,1:-1])
    u_out = np.sum(vel[-2,1:-1])
    #vel[0,1:-1] = vel[-1,1:-1]#*u_in/u_out
    vel[-1,1:-1] = vel[-2,1:-1]*u_in/u_out
    return

def walls_xmomentum(vel):
    vel[:,0] = - vel[:,1]
    vel[:,-1] = - vel[:,-2]
    return

def inflow_outflow_ymomentum(vel):
    vel[0,:] = vel[1,:]
    vel[-1,:] = vel[-2,:]
    return

def walls_ymomentum(vel):
    vel[:,0] = 0.0
    vel[:,-1] = 0.0
    return

def periodic_x(vel):
    vel[-1,1:-1] = vel[0,1:-1]
    return

def periodic_y(vel):
    vel[-1,:] = vel[1,:]
    vel[:,-1] = vel[:,0]
    return


inflow_outflow_xmomentum_noscale(u)
walls_xmomentum(u)
inflow_outflow_ymomentum(v)
walls_ymomentum(v)
#periodic_x(u)
#periodic_y(v)

time = 0
tic = timer.time()
start = timer.time()
tic = start
fig, (ax1, ax2) = plt.subplots(2,1, figsize=(6,12))
for steps in range(NUM_STEPS):
    
    #ut = np.roll(ut,1,0)
    #step 1, Calculate u_star and v_star
    
    J_u_x = 0.25*(u[:-1,1:-1] + u[1:,1:-1])**2
    J_u_y = 0.25*(u[1:-1,1:] + u[1:-1,:-1])*\
                 (v[2:-1,:] + v[1:-2,:]) 
                 
                
    J_u_x -= nu*(u[1:,1:-1]-u[:-1,1:-1])/dx
    J_u_y -= nu*(u[1:-1,1:]-u[1:-1,:-1])/dy
    
    ut[1:-1,1:-1] = u[1:-1,1:-1] - (DT/dxdy)*\
                    (J_u_x[1:,:]*dy - J_u_x[:-1,:]*dy + \
                     J_u_y[:,1:]*dx - J_u_y[:,:-1]*dx) + G*DT
    
    #np.roll(J_u_x,1,0) * dy
            
    J_v_x = 0.25*(u[:,2:-1]+u[:,1:-2])*(v[1:,1:-1]+v[:-1,1:-1])
    J_v_y = 0.25*(v[1:-1,1:]+v[1:-1,:-1])**2
    
    J_v_x -= nu*(v[1:,1:-1]-v[:-1,1:-1])/dx
    J_v_y -= nu*(v[1:-1,1:]-v[1:-1,:-1])/dy
    
    vt[1:-1,1:-1] = v[1:-1,1:-1] - (DT/dxdy)*\
                    (J_v_x[1:,:]*dy - J_v_x[:-1,:]*dy + \
                     J_v_y[:,1:]*dx - J_v_y[:,:-1]*dx)
    #vt = np.roll(vt,1,0)
            
    ##Boundaries intermediate
    inflow_outflow_xmomentum_noscale(ut)
    walls_xmomentum(ut)
    inflow_outflow_ymomentum(vt)
    walls_ymomentum(vt)
    #periodic_x(ut)
    #periodic_y(vt)
    
    

    #step 2, Calulate corrected presure field
    
    divergence = (ut[1:,1:-1]-ut[:-1,1:-1])/dx+\
                 (vt[1:-1,1:]-vt[1:-1,:-1])/dy

    prhs = RHO * divergence / DT
    for _ in range(50):
        p_next[:,:] = 0.0
        p_next[1:-1,1:-1] = (-prhs*dxdy**2 + \
                             (p[:-2,1:-1] + p[2:,1:-1])*dy**2 + 
                             (p[1:-1,:-2] + p[1:-1,2:])*dx**2)/\
                             (2.*(dx**2 + dy**2))
        p_next[0,:] = p_next[1,:] 
        p_next[-1,:] = p_next[-2,:]
        p_next[:,0] = p_next[:,1]
        p_next[:,-1] = p_next[:,-2]
        p[:,:] = p_next.copy()
    
    #step 3, Correct the velocity field (u and v)
    
    u[1:-1,1:-1] = ut[1:-1,1:-1] - (DT/RHO)*(p[2:-1,1:-1]-\
                                             p[1:-2,1:-1])/dx
    v[1:-1,1:-1] = vt[1:-1,1:-1] - (DT/RHO)*(p[1:-1,2:-1]-\
                                             p[1:-1,1:-2])/dy
    inflow_outflow_xmomentum_noscale(u)
    walls_xmomentum(u)
    inflow_outflow_ymomentum(v)
    walls_ymomentum(v)
    #periodic_x(u)
    #periodic_y(v)
    
    time += DT
    if (steps % PLOT_EVERY) == 0:
        toc = timer.time()
        divu = (u[1:,1:-1]-u[:1,1:-1])/dx + (v[1:-1,1:]-v[1:-1,:-1])/dy
        print(f'steps{steps} \t|divu|={np.linalg.norm(divu):.2e} \tTimes={time} \tSec/itr={(toc-tic)/PLOT_EVERY:.4f}s \tTotal Elapsed time={toc-start:.2f}s')
        uu = 0.5*(u[:-1,1:-1]+u[1:,1:-1])
        vv = 0.5*(v[1:-1,:-1]+v[1:-1,1:])
        ax1.clear()
        vel_fig = ax1.pcolormesh(xx,yy,np.sqrt(uu**2+vv**2), shading='auto', cmap='jet')
        #vel_fig = ax1.contourf(xx,yy, np.sqrt(uu**2+vv**2), levels=20, cmap='jet')
        vel_bar = plt.colorbar(vel_fig, ax=ax1)
        #ax1.quiver(xx,yy,uu,vv)
        ax1.set_xlim(0,LX)
        ax1.set_ylim(0,LY)
        ax1.set_aspect('equal')
        ax1.plot(
            int(NX/2)*dx + uu[int(NX/2),:],
            ynodes,
            'r-', lw=3
            )
        ax1.plot(
            int(NX/2)*dx+(0.1*ynodes*(LY-ynodes))/(2*MU),
            ynodes,
            'k-', lw=3
            )
        p_fig = ax2.pcolormesh(xx,yy, p[1:-1,1:-1], shading='auto', cmap='jet')
        p_bar = plt.colorbar(p_fig, ax=ax2)
        ax2.set_aspect('equal')
        plt.pause(0.01)
        vel_bar.remove()
        p_bar.remove()
        process = psutil.Process(os.getpid())
        print(f'Memory usage (RSS): {process.memory_info().rss/(1024**2):.2f} MB')
        tic = timer.time()
    plt.show()
