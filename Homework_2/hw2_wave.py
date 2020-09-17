#!/usr/bin/env python
#
# Homework 2
#
from __future__ import print_function
import os,sys

import numpy as np
print("")
print("numpy version: ",np.__version__)
print("")

import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 20, 5
plt.rcParams['lines.linewidth'] = 0.5

#import copy

# meshing parameters
dx = 0.1                    # grid size
x = np.arange(0,100.0+dx/2,dx)  # coordinates of grid points: 0-100
nx = len(x)                 # number of grid pints

# model parameters
rho = np.ones(nx)
kappa = np.ones(nx)

## material contrast
#x_discon = 60
#rho[int(x_discon/dx):nx] = 1.0     # rho contrast
#kappa[int(x_discon/dx):nx] = 4.0   # kappa contrast

c2 = kappa/rho                      # wavespeed square: c^2

# timing parameters
FACTOR = 0.25

dt = FACTOR * np.min(dx**2/c2)      # time step
dt = np.round(dt,decimals=4)        # round-off 4 digit

T_total = 10.0                      # simulation time: 0-200
t = np.arange(0,T_total+dt/2,dt)
nt = len(t)

print("time step dt = ",dt)
print("number of time steps nt = ",nt)
print("")

# initial condition
sigma = 0.1
# second-order system
u = np.exp(-sigma*(x-50)**2)        # displacement
v = np.zeros(u.shape)               # velocity 
u_old1 = np.copy(u)
u_new1 = np.copy(u)

# first-order system
T = kappa*(-2)*(x-50)*u*sigma       # stress
V = np.zeros(T.shape)               # velocity
T_old = np.copy(T)
V_old = np.copy(V)

# boundary contition = 'Dirichlet' or 'Neumann'
boundary = 'Neumann'

print("starting time loop with ... ",boundary)
print("")

# time marching
for it in np.arange(0,nt):

    ####################### 2nd order equation ################################
    #> TODO: implement your finite-difference scheme
    u_old1 = ..
    u = ..
    
    
    for i in np.arange(1,nx-1):
        u_new1[i] = ..

    # Dirichlet boundary condition
    if boundary == 'Dirichlet':
        ..

    # Neumann boundary condition
    if boundary == 'Neumann':
        ..

    #<TODO

    ####################### 1st order equations ###############################
    #> TODO: implement your finite-difference scheme
    V_old = ..
    for i in np.arange(1,nx-1):
        V[i] = ..
      
    T_old = ..
    for i in np.arange(1,nx-1):
        T[i] = ..
  
    # Dirichlet boundary condition
    if boundary == 'Dirichlet':
        ..

    # Neumann boundary condition
    if boundary == 'Neumann':
        ..

    ######################################################

    # figures
    if (it+1)%1000 == 0:
        time = it * dt
        time = np.round(time,decimals=4)        # round-off 4 digit
        print("time step ",it," - time = ",time)

        # plots velocity
        v = (u-u_old1)/dt  # velocity derived from displacement field
        plt.clf()
        plt.title('Time t = ' + str(time))
        plt.suptitle('Boundary condition = '+ boundary)
        plt.plot(x,V,color = 'red',label = '1st order')
        plt.plot(x,v,color = 'blue',label = '2nd order')
        plt.ylim(-1,1)
        plt.legend()
        #plt.show()

        # saves figure as pdf
        filename = 'figures/wave_plot_' + str(it) + '.pdf'
        plt.savefig(filename)
        print("plotted as ",filename)
        print("")

        # plots velocity difference
        #delta_v = V-v
        #plt.title('Velocity Differences - Time t = '+ str(time))
        #plt.plot(x,delta_v,label = 'velocity difference')
        #plt.ylim(-1,1)
        #plt.show()


print("")
print("done")
print("")
