#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 14:51:01 2025

@author: jasperwetzel
"""
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from git import Repo
import time as timer

Time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

plt.annotate(Time, xy=(0.7,0.95), xycoords='figure fraction', annotation_clip=False)

## Git repo

repo = Repo('.', search_parent_directories=True)
revsha = repo.head.object.hexsha[:8]
plt.annotate(f"[rev {revsha}]", xy=(0.05,0.95), xycoords='figure fraction', annotation_clip=False)


#values 
LY = 2e-3
LX = 50e-3
K = 0.4
CP = 3500
RHO = 1050
MU = 2.0e-3
nu = MU/RHO
U0 = 0.2
Q = 1.5e4
phi_0 = 1
phi_L = 0
T_LEFT = 298
T_initial = 298


#grid -Spacial 

NX = 5
NY = 5
dx = LX/NX
dy = LY/NY
x = np.linspace(0.5*dx, LX - 0.5*dx, NX)
y = np.linspace(0.5*dy, LY - 0.5*dy, NY)
xx, yy = np.meshgrid(x,y,indexing='ij')
T_max=[]
max_steps =[]
t_vals =[]  

#grid -temporal

sim_time = 0
steps = 0
t_END = 5
DT = 1e-8
PLOT_EVERY = 100
tic = timer.time()
start = timer.time()
tic = start

# velocity profile
def u_y(y):
    u_y = U0 * (3/2.) * (1-((2*y/LY)-1)**2)
    return u_y

ux = np.ones((NX+1,NY+2))
ux[0,1:-1] = ux[0,1:-1]*u_y(y)
ux[:,:] = np.tile(ux[0,:], (NX+1, 1))
ux[:,0] = 0.0
ux[:,-1] = 0.0

uy = np.ones((NX+2,NY+1))
uy[1:-1,0] = uy[1:-1,0]*u_y(y)
uy[:,:] = np.tile(uy[0,:], (NX+2, 1))
uy[0,:] = 0.0
uy[-1,:] = 0.0

#Temp Boundary conditions

def T_North(flux):
    flux[:,-1] = 0
    return

def f_South(flux):
    #temp[:,0] = -Q*2*dx/K - temp[:,1]
    flux[:,0] = -Q
    return

def T_West(temp,flux,flux_c):
    temp[0,1:-1] = T_LEFT*2 - temp[1,1:-1]
    flux[0,:] = 0.5*K * (temp[1,1:-1] - T_LEFT)/dx
    flux_c[0,:] = ux[-1,1:-1] * T_LEFT * RHO * CP
    return

def T_East(temp, flux):
    temp[-1,1:-1] = temp[-2,1:-1]
    flux[-1,:] = 0
    
    return

#Initalisation 

T = np.ones((NX+2,NY+2)) * T_initial

fig, ax1 = plt.subplots(figsize=(8,8))
#Temp_plot = plt.contour(T, origin ='lower', extent=(0,LX,0,LY),cmap ='hot')
#Temp_bar = plt.colorbar(Temp_plot, ax=ax1, label='Temprature (K)')

#fluxes 

T_EW = np.zeros((NX,NY))
T_NS = np.zeros((NX,NY))
T_EW_c = np.zeros((NX,NY))
T_NS_c = np.zeros((NX,NY))
x_flux = np.zeros((NX+1,NY))
y_flux = np.zeros((NX,NY+1))

#gs = fig.add_gridspec(1,1)
#ax0 = fig.add_subplot(gs[0,0])

T_North(T_EW)
f_South(T_EW)
T_North(T_NS)
f_South(T_NS)
f_South(T_NS_c)
T_East(T, T_EW)
T_West(T, T_EW, T_EW_c)

while sim_time < t_END:
    
    T_EW[:,:] = K * (T[1:,1:-1] - T[:-1,1:-1])/dx
    T_NS[:,:] = K * (T[1:-1,1:] - T[1:-1,:-1])/dy
    
    T_EW_c[1:-1,:] = 0.5*((T[2:-1,1:-1] + T[1:-2,1:-1]) * ux[1:-1,1:-1] -\
                  (T[1:-2,1:-1] + T[2:-1,1:-1]) * ux[1:-1,1:-1])*dy
        
    #T_NS_c[:,1:-1] = 0.5*((T[1:-1,2:-1] + T[1:-1,1:-2]) * v[1:-1,1:-1] -\
                  #(T[1:-1,1:-2] + T[1:-1,2:-1]) * v[1:-1,1:-1])*dx
    
    #Boundaries    
    #T_NS_c[:,0] = TSOUTH* v[1:-1,0] /dy
    #T_NS_c[:,-1] = TNORTH* v[1:-1,-1] /dy
    
    #T_EW_c[0,:] = u[0,1:-1] * (T[0,1:-1] + T[1,1:-1])*0.5 /dx
    #T_EW_c[-1,:] = u[-1,1:-1] * (T[-1,1:-1] + T[-2,1:-1])*0.5 /dx
    
    #T[1:-1,1:-1] = T[1:-1,1:-1] - DT*(((T_EW_c[1:,:]-T_EW_c[:-1,:]) +\
                                       #(T_NS_c[:,1:]-T_NS_c[:,:-1])) +\
                                      #((T_EW[1:,:]-T_EW[:-1,:]) +\
                                       #(T_NS[:,1:]-T_NS[:,:-1]))) /dx*dy
   
    #internal fluxes
    #x_flux[1:-1,:] = K * (T[2:-1, 1:-1] - T[1:-2, 1:-1]) / dx
    #K*(T[2:,1:-1]-2*T[1:-1,1:-1]-T[:-2,1:-1])/dx
                      #(0.5*((uy[1:-1,:-1]*T[2:,1:-1])-(uy[1:-1,1:]*T[:-2,1:-1]))*DT/dy)
                        
    #x_flux_D[1:-1,:] = (K*DT*(T[2:,1:-1]- 2*T[1:-1,1:-1]+T[:-2,1:-1])/dx**2)
    #x_flux_C[1:-1,:] = ux[:-1,1:-1]*(T[1:-2,1:-1])*RHO*CP
    
    #y_flux[:,1:-1] = K*(T[1:-1,2:-1]-T[1:-1,1:-2])/dy
                    #(0.5*((T[1:-1,2:]-T[1:-1,:-2]) * ux[:,1:-1]))*DT/dx
                    
    #y_flux_D[:,1:-1] = K*DT*(T[1:-1,2:]- 2*T[1:-1,1:-1]+T[1:-1,:-2])/dy
    #y_flux_C[:,1:-1] = 0#(0.5*(T[1:-1,2:]-T[1:-1,:-2]) * ux[:,1:-1])/dx
    
    #boundary 
    T_North(T_EW)
    f_South(T_EW)
    T_North(T_NS)
    f_South(T_NS)
    f_South(T_NS_c)
    f_South(T_NS_c)
    T_East(T, T_EW)
    T_West(T, T_EW, T_EW_c)
    
    T[1:-1,1:-1] = T[1:-1,1:-1] + (T_EW[:,:]+T_NS[:,:]) - (T_NS_c[:,:]+T_EW_c)
    #T[1:-1,1:-1] = T[1:-1,1:-1] + (x_flux[:,:] + y_flux[:,:])
    sim_time += DT
    steps +=1
    T_interior = T[1:-1, 1:-1]
    m = T_interior.max()
    T_max.append(m)
    max_steps.append(steps)

    t_vals.append(steps * DT)

    
    
    if (steps % PLOT_EVERY) == 0:
        toc = timer.time()
        divu = (T)/dx*dy
        print(f'steps{steps} \t|divu|={np.linalg.norm(divu):.2e} \tTimes={sim_time} \tSec/itr={(toc-tic)/PLOT_EVERY:.4f}s \tTotal Elapsed time={toc-start:.2f}s')
        ax1.clear()
        vel_fig = ax1.pcolormesh(xx,yy,T[1:-1,1:-1], shading='auto', cmap='hot')
        vel_bar = plt.colorbar(vel_fig, ax=ax1)
        ax1.set_xlim(0,LX)
        ax1.set_ylim(0,LY)
        ax1.set_title(f'Time: {sim_time:.3f} s')
        
        plt.pause(0.01)
        vel_bar.remove()
        
vel_bar = plt.colorbar(vel_fig, ax=ax1)     
plt.show()

fig, ax2 = plt.subplots(figsize=(7,4))
ax2.plot(t_vals, T_max, lw=2)
ax2.set_xlabel("Time [s]")
ax2.set_ylabel("Max temperature [K]")
ax2.set_title("Maximum temperature vs time")
ax2.grid(True)
plt.show()
        #ax0.cla()
        #im = ax0.imshow(T.T, origin= 'lower', extent=(0,LX,0,LY),cmap='hot')

        #plt.savefig(f'image_{n}')

        #temp_bar = plt.colorbar(im, ax=ax0, label='Temperature (K)')
        #plt.annotate(f"[rev {revsha}]", xy=(0.05,0.95), xycoords='figure fraction', annotation_clip=False)

        #plt.annotate(Time, xy=(0.7,0.95), xycoords='figure fraction', annotation_clip=False)

        #plt.grid()
        #ax0.clear()
        #plt.pause(0.01)
        #temp_bar.remove()


#plt.plot(g, T_HX, 'bo')
#plt.xlabel("Position (m)")
#plt.ylabel("Temprature (K)")
#plt.title("Profile lines of temperature along Hx/2")

#fig_2 = plt.figure(figsize=(8,8))
#plt.grid()
#ax0.clear()

#plt.plot(g, T_HY, 'ro')
#plt.xlabel("Position (m)")
#plt.ylabel("Temprature (K)")
#plt.title("Profile lines of temperature along Hy/2")


