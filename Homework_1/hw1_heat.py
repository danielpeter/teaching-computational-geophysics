#!/usr/bin/env python
#
# Homework 1 - Heat diffusion
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
dx = 1.0                    # grid size
L  = 100.0                  # length
x = np.arange(0,L+dx/2,dx)  # coordinates of grid points: 0-100
nx = len(x)                 # number of grid points

# model parameters (diffusivity)
D = 1.0

# timing parameters
FACTOR = 0.3

print("heat diffusion - finite-difference scheme")
print("  factor = ",FACTOR)

# scales time step size
dt = FACTOR * dx**2 / D             # time step
dt = np.round(dt,decimals=4)        # round-off 4 digit

Time = 25.0                         # simulation time
t = np.arange(0,Time+dt/2,dt)
nt = len(t)

print("")
print("  time step dt = ",dt)
print("  space     dx = ",dx)
print("")
print("  number of time steps nt = ",nt)
print("  number of grid points   = ",nx)
print("")

# initial condition (explicit loop)
Te = np.zeros(nx)
#Te[50] = 1.0

# initial condition (implicit loop)
Ti = np.zeros(nx)
#Ti[50] = 1.0

# matrix for the system of linear equations: A T_n = T_n-1 for the implicit scheme
A = np.zeros((nx,nx))
# Coefficients b_1 & b_N
A[0,0] = 1.0
A[nx-1,nx-1] = 1.0
for i in np.arange(1,nx-1):
    #> TODO: implement your finite-difference scheme
    ..               # a_i
    ..               # b_i
    ..               # c_i
    #<TODO

# analytical solution
sigma = 1.0
K = 1.0
Tmax = 1.0
T_analytical = np.zeros(nx)

# initial temperature profile (analytical solution)
time = 0.0
T0 = np.zeros(nx)
for i in np.arange(0,nx):
    T0[i] = Tmax / np.sqrt(1.0 + 4.0 * time * K/sigma**2) * np.exp( -(x[i] - 50.0)**2 /(sigma**2 + 4.0 * time * K) )

# sets as initial Temperature profile
Te[:] = np.copy(T0)
Ti[:] = np.copy(T0)

# time marching
for it in np.arange(0,nt):

    ####################### explicit scheme ################################
    #> TODO: implement your finite-difference scheme

    # update the temperature field
    ..

    # formula for explicit scheme
    ..

    # boundary conditions
    ..

    #<TODO

    ####################### implicit scheme ###############################
    #> TODO: implement your finite-difference scheme

    # update the temperature field
    ..

    # solve linear system: A Ti = Ti_old  for unknowns in Ti
    ..

    #<TODO

    ######################################################

    # figures
    if it % 10 == 0:
        time = (it+1) * dt
        time = np.round(time,decimals=4)        # round-off 4 digit
        print("time step ",it," out of ",nt," - time = ",time)

        # analytical solution
        for i in np.arange(0,nx):
            T_analytical[i] = Tmax / np.sqrt(1.0 + 4.0 * time * K / sigma**2) * np.exp( -(x[i] - 50.0)**2 / (sigma**2 + 4.0 * time * K) )

        # plots velocity
        plt.clf()
        fig = plt.figure()
        plt.title("Time t = " + str(time))
        plt.text(x=0.8, y=0.9, s="dt = {} * dx**2 / D".format(FACTOR), fontsize=12, ha="center", transform=fig.transFigure)
        plt.plot(x,Te,color = 'red',label = 'explicit scheme')
        plt.plot(x,Ti,color = 'blue',label = 'implicit scheme')
        plt.plot(x,T_analytical,color = 'black',label = 'analytical solution')
        plt.ylim(0,1)
        plt.legend()

        # saves figure as pdf
        filename = 'figures/heat_plot_' + str(it) + '.pdf'
        plt.savefig(filename)
        #plt.show()

        print("plotted as ",filename)
        print("")

print("")
print("done")
print("")
