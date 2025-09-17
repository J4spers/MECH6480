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

step_vel_0 = np.genfromtxt(f'exper-xvel-re230-xh0.xy', names=True, delimiter=' ')
step_vel_6 = np.genfromtxt(f'exper-xvel-re230-xh6.xy', names=True, delimiter=' ')

#paramiters
LX = 0.6    #m
LY = 0.03   #m

RHO = 998.0     #kg/m^3
MU = 0.001      #Pa.s
nu = MU/RHO
RE = 230

#Boundary conditions

STEP_HEIGHT = 0.01
STEP_LENGTH = 0.22
PLOT_LOC = [0.22,0.28]
u_inlet = RE*nu/STEP_HEIGHT

# grid - spacial

NY = 30
NX = int(NY*(LX/LY))
dx = LX/NX
dy = LY/NY
dxdy = dx*dy

step_height = int(STEP_HEIGHT/dy)
step_width  = int(STEP_LENGTH/dx)
print(step_height,step_width)
PLOT_IDX = [int(PLOT_LOC[0]/dx),int(PLOT_LOC[1]/dx)]

xnodes = np.linspace(dx/2., LX-dx/2., NX)
ynodes = np.linspace(dy/2., LY-dy/2., NY)
xx, yy = np.meshgrid(xnodes, ynodes, indexing='ij')

#grid - temporal 
DT = 0.005
#TOTAL_TIME = 100
#num_steps = int(TOTAL_TIME/DT)
NUM_STEPS = 5000
PLOT_EVERY = 100
N_PRESSURE_ITS = 500

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
p_next = np.zeros_like(p)

uu = np.zeros((NX,NY))
vv = np.zeros((NX,NY))

#Initialisation 
u[:,step_height+1:] = u_inlet

#Boundary conditions

def inflow_outflow_xmomentum(vel):
    vel[0,(step_height+1):-1] = u_inlet
    u_in =np.sum(vel[0,(step_height+1):-1])
    u_out = np.sum(vel[-2,1:-1])
    vel[-1,1:-1] = vel[-2,1:-1]#*u_in/u_out
    return

def walls_xmomentum(vel):
    vel[(step_width+1):,0] = - vel[(step_width+1),1]
    vel[:,-1] = - vel[:,-2]
    
    vel[:(step_width+1),:(step_height+1)] = 0.0 # internal and right wall = 0
    vel[:(step_width+1),(step_height)] = -vel[:(step_width+1),(step_height+1)] #top of step
    return

def inflow_outflow_ymomentum(vel):
    vel[0,:] = vel[1,:]
    vel[-1,:] = vel[-2,:]
    return

def walls_ymomentum(vel):
    vel[(step_width+1):,0] = 0.0 #bottom wall
    vel[:,-1] = 0.0             #top wall
    
    vel[:(step_width+1),:(step_height+1)] = 0.0 # inside and top of step
    vel[step_width,:(step_height+1)] = -vel[(step_width+1), :(step_height+1)] #front of step
    return

inflow_outflow_xmomentum(u)
walls_xmomentum(u)
inflow_outflow_ymomentum(v)
walls_ymomentum(v)

time = 0
tic = timer.time()
start = timer.time()
tic = start
fig, (ax1, ax2) = plt.subplots(2,1, figsize=(6,12))
for steps in range(NUM_STEPS):
    
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
    inflow_outflow_xmomentum(ut)
    walls_xmomentum(ut)
    inflow_outflow_ymomentum(vt)
    walls_ymomentum(vt)
    

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
        p_next[0,(step_height+1):] = p_next[1,(step_height+1):] 
        p_next[-1,:] = -p_next[-2,:]
        p_next[(step_width+1):,0] = p_next[(step_width+1):,1] #bottom wall after step
        p_next[:(step_width+1),:(step_height+1)] = 0.0 # pressure inside step to zero
        p_next[:(step_width+1),step_height] = p_next[:(step_width+1),(step_height+1)] # top of step
        p_next[step_width, :(step_height+1)] = p_next[(step_width+1), :(step_height+1)]
        p_next[:,-1] = p_next[:,-2]
        p[:,:] = p_next.copy()
    
    #step 3, Correct the velocity field (u and v)
    
    u[1:-1,1:-1] = ut[1:-1,1:-1] - (DT/RHO)*(p[2:-1,1:-1]-\
                                             p[1:-2,1:-1])/dx
    v[1:-1,1:-1] = vt[1:-1,1:-1] - (DT/RHO)*(p[1:-1,2:-1]-\
                                             p[1:-1,1:-2])/dy
    inflow_outflow_xmomentum(u)
    walls_xmomentum(u)
    inflow_outflow_ymomentum(v)
    walls_ymomentum(v)
    
    time += DT
    if (steps % PLOT_EVERY) == 0:
        toc = timer.time()
        divu = (u[1:,1:-1]-u[:1,1:-1])/dx + (v[1:-1,1:]-v[1:-1,:-1])/dy
        print(f'steps{steps} \t|divu|={np.linalg.norm(divu):.2e} \tTimes={time} \tSec/itr={(toc-tic)/PLOT_EVERY:.4f}s \tTotal Elapsed time={toc-start:.2f}s')
        uu = 0.5*(u[:-1,1:-1]+u[1:,1:-1])
        vv = 0.5*(v[1:-1,:-1]+v[1:-1,1:])
        ax1.clear()
        umag = np.sqrt(uu**2+vv**2)
        #umag[]
        vel_fig = ax1.pcolormesh(xx,yy,umag, shading='auto', cmap='jet')
        #vel_fig = ax1.contourf(xx,yy, np.sqrt(uu**2+vv**2), levels=20, cmap='jet')
        vel_bar = plt.colorbar(vel_fig, ax=ax1)
        #ax1.quiver(xx,yy,uu,vv)
        ax1.set_xlim(0,LX)
        ax1.set_ylim(0,LY)
        ax1.set_aspect('equal')
        ax1.plot(
            PLOT_LOC[0] + uu[PLOT_IDX[0],:],
            ynodes,
            'k-', lw=3
            )
        ax1.plot(
            PLOT_LOC[1] + uu[PLOT_IDX[1],:],
            ynodes,
            'k-', lw=3
            )
        ax1.fill_between(
            [0, STEP_LENGTH],
            0,
            STEP_HEIGHT,
            color='white'
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

uu = 0.5*(u[:-1,1:-1]+u[1:,1:-1])
vv = 0.5*(v[1:-1,:-1]+v[1:-1,1:])
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(12,6))
ax1.plot(
    step_vel_0['y'],
    step_vel_0['ux'], 'ro', label='experimental')
ax1.plot(
    ynodes,
    uu[PLOT_IDX[0],:], 'b-', label='simulation')
ax1.setxlabel('y (m)')
ax1.setylabel('u (m/s)')
ax2.plot(
    step_vel_6['y'],
    step_vel_6['ux'], 'ro', label='experimental')
ax2.plot(
    ynodes,
    uu[PLOT_IDX[1],:], 'b-', label='simulation')
ax2.setxlabel('y (m)')
ax2.setylabel('u (m/s)')

