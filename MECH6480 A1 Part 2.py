#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 13:16:01 2025

@author: jasperwetzel
"""

import numpy as np
import matplotlib.pyplot as plt
import time as timer
import psutil, os
from datetime import datetime
from git import Repo
plt.close('all')

#udata = np.genfromtxt('exper-xvel-re230-xh0.xy', names=True, delimiter=' ')
#vdata = np.genfromtxt('exper-xvel-re230-xh6.xy', names=True, delimiter=' ')

Time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

## Git repo

repo = Repo('.', search_parent_directories=True)
revsha = repo.head.object.hexsha[:8]

#plt.annotate(f"[rev {revsha}]", xy=(0.05,0.95), xycoords='figure fraction', annotation_clip=False)
#plt.annotate(Time, xy=(0.7,0.95), xycoords='figure fraction', annotation_clip=False)

#paramiters
LX = 2*np.pi
LY = 2*np.pi

RHO = 1
nu = 0.001
MU = nu*RHO

# grid - spacial

NX = 16 #5
NY = 16 #5
dx = LX/NX
dy = LY/NY
dxdy = dx*dy

xnodes = np.linspace(dx/2., LX-dx/2., NX)
ynodes = np.linspace(dy/2., LY-dy/2., NY)
xx, yy = np.meshgrid(xnodes, ynodes, indexing='ij')

ux = np.linspace(0,LX,NX+1) 
vy = np.linspace(0,LY,NY+1) 

uxx, uyy = np.meshgrid(ux,ynodes, indexing='ij') # positions of u nodes in the x and y
vxx, vyy = np.meshgrid(xnodes, vy, indexing='ij') # position of v nodes in the x and y

#grid - temporal 

DT = 0.01
TOTAL_TIME = 200
num_steps = int(TOTAL_TIME/DT)
#NUM_STEPS = 5000
PLOT_EVERY = 1000
N_PRESSURE_ITS = 20

def UIN(x,y):
    uin = np.sin(x) * np.cos(y)
    return uin
def VIN(x,y):
    vin = - np.cos(x) * np.sin(y)
    return vin

#ana sol
def u_sol(x,y,t):
    u_sol = np.sin(x) * np.cos(y) * np.exp(-2*nu*t)
    return u_sol
def v_sol(x,y,t):
    v_sol = -1* np.cos(x) * np.sin(y) * np.exp(-2*nu*t)
    return v_sol

def order_of_accuracy(e_course, e_fine, r):
    return np.log(e_course/e_fine)/np.log(r)

#alocate space and initilisation

u = np.ones((NX+1, NY+2))
u[:,1:-1] = UIN(uxx,uyy)

uA = np.ones_like(u)
uA[:,1:-1] = UIN(uxx,uyy)
ut = np.zeros_like(u)
J_u_x = np.zeros((NX,NY))
J_u_y = np.zeros((NX,NY+1))

v = np.ones((NX+2, NY+1))
v[1:-1,:] = VIN(vxx,vyy)
vA = np.ones_like(v)
vA[1:-1,:] = VIN(vxx,vyy)
vt = np.zeros_like(v)
J_v_x = np.zeros((NX+1,NY))
J_v_y = np.zeros((NX,NY))

p = np.zeros((NX+2,NY+2))
prhs = np.zeros_like(p)
p_next = np.zeros_like(p)

uu = np.zeros((NX,NY))
vv = np.zeros((NX,NY))

#Boundary conditions

def u_update_east(vel):
    vel[:,-1] = vel[:,1]
    vel[:,0] = vel[:,-2]
    return

def u_update_north(vel):
    vel[0,:] = vel[-1,:]
    return

def v_update_west(vel):
    vel[-1,:] = vel[1,:]
    vel[0,:] = vel[-2,:]
    
    return

def v_update_south(vel):
    vel[:,-1] = vel[:,0]
    return

u_update_east(u)
u_update_north(u)
v_update_west(v)
v_update_south(v)

time = 0
tic = timer.time()
start = timer.time()
tic = start
fig, (ax1, ax2) = plt.subplots(2,1, figsize=(6,12))
fig2, ax3 = plt.subplots(figsize=(6,6))
for steps in range(num_steps):
    
    #step 1, Calculate u_star and v_star
    
    J_u_x = 0.25*(u[:-1,1:-1] + u[1:,1:-1])**2
    J_u_y[:,:] = 0.25*(u[:-1,1:] + u[:-1,:-1])*\
                 (v[1:-1,:] + v[:-2,:]) 

    J_u_x -= nu*(u[1:,1:-1]-u[:-1,1:-1])/dx
    J_u_y -= nu*(u[:-1,1:] - u[:-1,:-1])/dy
    
    ut[1:-1,1:-1] = u[1:-1,1:-1] - (DT/dxdy)*\
                    (J_u_x[1:,:]*dy - J_u_x[:-1,:]*dy + \
                     J_u_y[1:,1:]*dx - J_u_y[1:,:-1]*dx)
    
    ut[0,1:-1] = u[0,1:-1] - (DT/dxdy)*\
                    (J_u_x[0,:]*dy - J_u_x[-1,:]*dy + \
                     J_u_y[0,1:]*dx - J_u_y[0,:-1]*dx)
            
            
    J_v_x[:,:] = 0.25*(u[:,1:-1]+u[:,:-2])*\
                 (v[1:,:-1]+v[:-1,:-1])
    
    J_v_y = 0.25*(v[1:-1,1:]+v[1:-1,:-1])**2
    
    
    J_v_x -= nu*(v[1:,:-1]-v[:-1,:-1])/dx
    
    J_v_y -= nu*(v[1:-1,1:]-v[1:-1,:-1])/dy
    
    
    vt[1:-1,1:-1] = v[1:-1,1:-1] - (DT/dxdy)*\
                    (J_v_x[1:,1:]*dy - J_v_x[:-1,1:]*dy + \
                     J_v_y[:,1:]*dx - J_v_y[:,:-1]*dx)
    
    vt[1:-1,0] = v[1:-1,0] - (DT/dxdy)*\
                    (J_v_x[1:,0]*dy - J_v_x[:-1,0]*dy + \
                     J_v_y[:,0]*dx - J_v_y[:,-1]*dx)
    
    
    ##Boundaries intermediate
    u_update_east(ut)
    u_update_north(ut)
    v_update_west(vt)
    v_update_south(vt)
    

    #step 2, Calulate corrected presure field
    
    divergence = (ut[1:,1:-1]-ut[:-1,1:-1])/dx+\
                 (vt[1:-1,1:]-vt[1:-1,:-1])/dy

    prhs = RHO * divergence / DT
    prhs -= prhs.mean()
    for _ in range(N_PRESSURE_ITS):
        p_next[:,:] = 0.0
        p_next[1:-1,1:-1] = (-prhs*dxdy**2 + \
                             (p[:-2,1:-1] + p[2:,1:-1])*dy**2 + 
                             (p[1:-1,:-2] + p[1:-1,2:])*dx**2)/\
                             (2.*(dx**2 + dy**2))
        p_next[0,:] = p_next[-2,:] 
        p_next[-1,:] = p_next[1,:]
        p_next[:,0] = p_next[:,-2]
        p_next[:,-1] = p_next[:,1]
        p[:,:] = p_next.copy()
    
    #step 3, Correct the velocity field (u and v)
    
    #u[:,1:-1] = ut[:,1:-1] - (DT/RHO)*(p[1:,1:-1]-\
                                             #p[:-1,1:-1])/dx
        
    #v[1:-1,:] = vt[1:-1,:] - (DT/RHO)*(p[1:-1,1:]-\
                                             #p[1:-1,:-1])/dy

        
    u[:-1,1:-1] = ut[:-1,1:-1] - (DT/RHO)*(p[1:-1,1:-1]-\
                                             p[:-2,1:-1])/dx
        
    v[1:-1,:-1] = vt[1:-1,:-1] - (DT/RHO)*(p[1:-1,1:-1]-\
                                             p[1:-1,:-2])/dy
    u_update_east(u)
    u_update_north(u)
    v_update_west(v)
    v_update_south(v)
    
    time += DT
    #Anylitical solution
    uA[:,1:-1] = u_sol(uxx, uyy, time)
    vA[1:-1,:] = v_sol(vxx, vyy, time)

    if (steps % PLOT_EVERY) == 0:
        toc = timer.time()
        divu = (u[1:,1:-1]-u[:-1,1:-1])/dx + (v[1:-1,1:]-v[1:-1,:-1])/dy
        print(f'steps{steps} \t|divu|={np.linalg.norm(divu):.2e} \tTimes={time} \tSec/itr={(toc-tic)/PLOT_EVERY:.4f}s \tTotal Elapsed time={toc-start:.2f}s')
        uu = 0.5*(u[:-1,1:-1]+u[1:,1:-1]) #interpolate to cell centred velocities for u
        vv = 0.5*(v[1:-1,:-1]+v[1:-1,1:]) #interpolate to cell centred velocities for v
        ax1.clear()
        vel_fig = ax1.pcolormesh(xx,yy,np.sqrt(uu**2+vv**2), shading='auto', cmap='jet')
        #vel_fig = ax1.contourf(xx,yy, np.sqrt(uu**2+vv**2), levels=20, cmap='jet')
        vel_bar = plt.colorbar(vel_fig, ax=ax1)
        ax1.quiver(xx,yy,uu,vv)
        ax1.set_xlim(0,LX)
        ax1.set_ylim(0,LY)
        ax1.set_aspect('equal')
        plt.annotate(f"[rev {revsha}]", xy=(0.05,0.95), xycoords='figure fraction', annotation_clip=False)
        plt.annotate(Time, xy=(0.7,0.95), xycoords='figure fraction', annotation_clip=False)
        uAA = 0.5*(uA[:-1,1:-1]+uA[1:,1:-1]) #cell centred anylitical velocities for u 
        vAA = 0.5*(vA[1:-1,:-1]+vA[1:-1,1:]) #cell centred anylitical velocities for v
        
        ax1.plot(
            np.pi + np.sqrt(uu[int(np.pi/dx),:]**2+vv[int(np.pi/dx),:]**2),
            ynodes,
            'r-', lw=3
            )
        #anaylitical sol
        ax1.plot(
            np.pi + np.sqrt(uAA[int(np.pi/dx),:]**2+vAA[int(np.pi/dx),:]**2),
            ynodes,
            'k-', lw=3
            )
        # changed pressure graph to anylitical vel
        
        real_sol_fig = ax2.pcolormesh(xx,yy, np.sqrt(uAA**2+vAA**2), shading='auto', cmap='jet')
        p_fig = ax2.pcolormesh(xx,yy, np.sqrt(uAA**2+vAA**2), shading='auto', cmap='jet')
        p_bar = plt.colorbar(p_fig, ax=ax2)
        ax2.set_xlabel("x [m]")
        ax2.set_ylabel("y [m]")
        
        ax1.set_ylabel("y [m]")
        ax1.set_title("Calculated velocity spirals")
        ax2.set_title("Analytical velocity spirals")
        ax2.quiver(xx,yy,uAA,vAA)
        ax2.set_aspect('equal')
        
        
        plt.pause(0.01)
        vel_bar.remove()
        p_bar.remove()
        
        
        process = psutil.Process(os.getpid())
        print(f'Memory usage (RSS): {process.memory_info().rss/(1024**2):.2f} MB')
        tic = timer.time()
        
        
        vel_dif_fig = ax3.pcolormesh(xx,yy, (np.sqrt(uu**2+vv**2)-np.sqrt(uAA**2+vAA**2)), shading='auto', cmap='jet')
        ax3.set_xlabel("x [m]")
        ax3.set_ylabel("y [m]")
        ax3.set_title("Velocity Magnitude, Analytical - Calculated")
        plt.annotate(f't={time:.3f}s', xy=(0.4,0.95), xycoords='figure fraction', annotation_clip=False)

        dif_bar = plt.colorbar(vel_dif_fig, ax=ax3) 
        plt.pause(0.1)
        dif_bar.remove()
        ax3.clear()
    plt.show()
dif_bar = plt.colorbar(vel_dif_fig, ax=ax3)
p_bar = plt.colorbar(p_fig, ax=ax2)
vel_bar = plt.colorbar(vel_fig, ax=ax1)
plt.annotate(f"[rev {revsha}]", xy=(0.06,0.05), xycoords='figure fraction', annotation_clip=False)
plt.annotate(Time, xy=(0.7,0.95), xycoords='figure fraction', annotation_clip=False)
#plt.annotate(f't={time:.3f}s', xy=(0.4,0.85), xycoords='figure fraction', annotation_clip=False)
error = np.sqrt(uAA[int(np.pi/dx),:]**2+vAA[int(np.pi/dx),:]**2) - np.sqrt(uu[int(np.pi/dx),:]**2+vv[int(np.pi/dx),:]**2)
plt.annotate(f't={time:.3f}s', xy=(0.4,0.95), xycoords='figure fraction', annotation_clip=False)


N12 = [ 0.02568792,  0.01054522, -0.01530462, -0.03421679, -0.04116228,
       -0.03327622,  0.03477301,  0.04403187,  0.02524241, -0.00497241,
       -0.02987671, -0.03046092]

N16 = []

N32 = [-0.51713241, -0.19901926, -0.05983263,  0.05947807,  0.14434592,
        0.19406395,  0.21409115,  0.21463346,  0.20271344,  0.18335199,
        0.15844606,  0.13009066,  0.09484183,  0.04687136, -0.04509656,
       -0.19938757, -0.23949757, -0.16322735, -0.10460633, -0.04514011,
        0.01066744,  0.06070767,  0.10360805,  0.13798922,  0.16367511,
        0.17770964,  0.17490703,  0.14667171,  0.07941861, -0.04475341,
       -0.23133964, -0.57291811]

N64 = [-0.96970402, -0.77740722, -0.71091142, -0.5869667 , -0.45315975,
       -0.31165217, -0.18004369, -0.06195945,  0.04418897,  0.14058225,
        0.21403281,  0.27300009,  0.30856791,  0.33091039,  0.33202309,
        0.32519702,  0.32489792,  0.29061037,  0.21793011,  0.18660191,
        0.1796937 ,  0.16384306,  0.07795128, -0.14816365, -0.38485614,
       -0.50186798, -0.50673922, -0.44755123, -0.36044551, -0.29290606,
       -0.23545012, -0.18203913, -0.13962217, -0.13604675, -0.16678306,
       -0.21701432, -0.29752975, -0.36974868, -0.42389169, -0.43414914,
       -0.35004039, -0.13587041,  0.17260596,  0.41371648,  0.44738686,
        0.40844794,  0.38764619,  0.37898769,  0.37709343,  0.35311317,
        0.32283907,  0.27950574,  0.22201819,  0.15488364,  0.07087871,
       -0.01732185, -0.12082343, -0.22342564, -0.33021151, -0.44291723,
       -0.54549912, -0.6525621 , -0.71415573, -0.71602205]

n12 = sum(N12)/len(N12)
n16 = sum(N16)/len(N16)
n32 = sum(N32)/len(N32)
n64 = sum(N64)/len(N64)

print(order_of_accuracy(n12,n16,(16/12)))