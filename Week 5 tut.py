#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 15:05:48 2025

@author: jasperwetzel
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from git import Repo
from matplotlib . animation import FuncAnimation

Time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

plt.annotate(Time, xy=(0.7,0.95), xycoords='figure fraction', annotation_clip=False)

## Git repo
repo = Repo('.', search_parent_directories=True)
revsha = repo.head.object.hexsha[:8]
plt.annotate(f"[rev {revsha}]", xy=(0.05,0.95), xycoords='figure fraction', annotation_clip=False)


#values 

LY = 10
LX = 10
ALPHA = 4


#grid -Spacial 

NX = 21
NY = 21
dx = LX/NX
dy = LY/NY
x = np.linspace(0.5*dx, LX - 0.5*dx, NX)
y = np.linspace(0.5*dy, LY - 0.5*dy, NY)
xx, yy = np.meshgrid(x,y,indexing='ij')

#grid -temporal
n = -1
sim_time = 0
iteration =0
t_END = 0.5
DT = 0.001
PLOT_EVERY = 500

# initial conditions and boundary conditions 
T_LEFT = 400
T_TOP = 350
T_initial = 300

T = np.ones((NX,NY)) * T_initial
T[0,:-1] = T_LEFT
T[1:,-1] = T_TOP


plt.imshow(T.T, origin ='lower', extent=(0,LX,0,LY),cmap ='hot')
plt.colorbar(label='Temprature (K)')


#fluxes 

x_flux = np.zeros((NX+1,NY))
y_flux = np.zeros((NX,NY+1))

fig = plt.figure(figsize=(8,8))
gs = fig.add_gridspec(1,1)
ax0 = fig.add_subplot(gs[0,0])

while sim_time < t_END:
    #internal fluxes
    x_flux[1:-1,:] = ALPHA * (T[1:,:]-T[:-1,:])/dx
    y_flux[:,1:-1] = ALPHA * (T[:,1:]-T[:,:-1])/dy
    
    
    #boundary
    x_flux[-1,:] = 0
    y_flux[:,0] = 0
    x_flux[0,:] = ALPHA * (T[1,:] - T_LEFT)/(dx/2.)
    y_flux[:,-1] = ALPHA * (T[:,-2] - T_TOP)/(dy/2.)
    T[0,:-1] = T_LEFT
    T[1:,-1] = T_TOP
    
    T = T + DT*(dy*x_flux[1:,:] - dy*x_flux[:-1,:] + dx*y_flux[:,1:] - dx*y_flux[:,:-1])/(dx*dy)
    sim_time += DT
    iteration +=1
    
    if np.mod(iteration, PLOT_EVERY) == 0 :
        n += 1
        ax0.cla()
        im = ax0.imshow(T.T, origin= 'lower', extent=(0,LX,0,LY),cmap='hot')
        ax0.set_title(f'Time: {sim_time:.3f} s')
        plt.savefig(f'image_{n}')

plt.colorbar(im, ax=ax0, label='Temperature (K)')
plt.annotate(f"[rev {revsha}]", xy=(0.05,0.95), xycoords='figure fraction', annotation_clip=False)

plt.annotate(Time, xy=(0.7,0.95), xycoords='figure fraction', annotation_clip=False)

nframes = 99
plt. subplots_adjust (top=1, bottom =0, left =0, right =1)

def animate(i):
    im = plt.imread(f'image_{i}.png')
    plt.imshow(im)
anim = FuncAnimation (plt.gcf(), animate , frames=nframes)
anim.save('output.gif', writer='pillow')


fig_1 = plt.figure(figsize=(8,8))
plt.grid()
ax0.clear()
T_HX = (T[:,10] + T[:,11])/2
T_HY = (T[10,:] + T[11,:])/2
g = np.linspace(0, LX, num=NX)

plt.plot(g, T_HX, 'bo')
plt.xlabel("Position (m)")
plt.ylabel("Temprature (K)")
plt.title("Profile lines of temperature along Hx/2")

fig_2 = plt.figure(figsize=(8,8))
plt.grid()
ax0.clear()
T_HX = (T[:,10] + T[:,11])/2
T_HY = (T[10,:] + T[11,:])/2

plt.plot(g, T_HY, 'ro')
plt.xlabel("Position (m)")
plt.ylabel("Temprature (K)")
plt.title("Profile lines of temperature along Hy/2")

plt.show()