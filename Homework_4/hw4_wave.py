#!/usr/bin/env python
#
# Homework 4
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

# meshing parameters
dx = 0.1                    # grid size
L = 100.0                   # length
x = np.arange(0,L,dx)       # coordinates of grid points: 0-100

# makes sure nx is even
if len(x)%2 != 0: x = np.arange(0,L+dx,dx)
#power of 2 length 1024, better for FFT
#x = np.arange(0,102.4+dx,dx)

# number of grid points
nx = len(x)

# wavenumber increment
dk = 2.0 * np.pi / nx / dx

print("")
print("grid step dx = ",dx);
print("number of grid steps nx = ",nx)
print("")

# model parameters
rho = np.ones(nx)
kappa = np.ones(nx)

x_discon = 60
rho[int(x_discon/dx):nx] = 1.0     # rho contrast
kappa[int(x_discon/dx):nx] = 4.0   # kappa contrast

# wavespeed square: c^2
c2 = kappa / rho            

# timing parameters
## Pseudo-spectral (PS) time step size
# TODO: use different factors to find good stability
dt = 0.1 * min(dx / np.sqrt(c2) / (2.0 * np.pi))

dt = np.round(dt,6)
print("PS time step dt = ",dt)

## Finite-difference (FD) time step (see homework 2)
do_compare_with_FD = True
if do_compare_with_FD:
    FD_dt = 0.3 * min(dx**2 / c2)
    FD_dt = np.round(FD_dt,6)
    print("FD time step dt = ",FD_dt)
    # using smaller time step as common for both schemes
    dt = min(dt,FD_dt)
    FD_dt = dt
print("")

## simulation time: 50 sec
t = np.arange(0,50,dt)

nt = len(t)

print("using time step dt = ",dt)
print("number of time steps nt = ",nt)
print("")

# initial condition
sigma = 0.1
u = np.exp(- sigma * (x-50)**2)        # displacement

# for (velocity,stress) system:
# instead of displacement, we can set stress (derived from displacement)
T = kappa * (-sigma * 2 * (x-50)) * u  # stress: T = kappa * du_dx
V = np.zeros(len(T))                   # velocity

V_old1 = np.copy(V)
V_old2 = np.copy(V)

T_old1 = np.copy(T)
T_old2 = np.copy(T)

## FD scheme
if do_compare_with_FD:
    FD_T = np.copy(T)
    FD_T_old = np.copy(FD_T)
    FD_V = np.copy(V)
    FD_V_old = np.copy(FD_V)

# time marching
for it in np.arange(0,nt):
    # current time
    t = dt * it

    ## Pseudo-spectral scheme
    # velocity
    V_old2 = np.copy(V_old1)
    V_old1 = np.copy(V)

    # stress
    T_old2 = np.copy(T_old1)
    T_old1 = np.copy(T)

    # displacement field (for comparisons/plotting only)
    u = u + dt * V;

    #>TODO: implement your pseudo-spectral method
    #       based on a Fourier transform (fft)
    #       (note: if X is real, then Y = fft(X) is conjugate symmetric
    #        and the number of unique points is ceil((n+1)/2)
    ..

    # Dirichlet boundary condition
    ..

    # Neumann boundary condition
    ..

    #<TODO

    ## FD scheme (see homework 2)
    if do_compare_with_FD:
        FD_V_old = np.copy(FD_V)
        FD_T_old = np.copy(FD_T)
        FD_V[1:nx-1] = FD_V_old[1:nx-1] + FD_dt / 2 / dx / rho[1:nx-1]   * (FD_T_old[2:nx]-FD_T_old[0:nx-2])
        FD_T[1:nx-1] = FD_T_old[1:nx-1] + FD_dt / 2 / dx * kappa[1:nx-1] * (FD_V_old[2:nx]-FD_V_old[0:nx-2])
        # Dirichlet boundary condition
        FD_V[0] = 0.0             # left
        FD_T[0] = FD_T[1]
        #FD_V[nx-1] = 0.0         # right
        #FD_T[nx-1] = FD_T[nx-2]
        # Neumann boundary condition
        #FD_V[0] = FD_V[1]        # left
        #FD_T[0] = 0.0
        FD_V[nx-1] = FD_V[nx-2]   # right
        FD_T[nx-1] = 0.0

    # figures
    if (it+1)%int(nt/10) == 0:
        time = it * dt
        time = np.round(time,decimals=4)        # round-off 4 digit
        print("time step ",it," - time = ",time)

        # plots velocity
        plt.clf()
        plt.title('Time t = ' + str(time))
        plt.plot(x,V,color = 'red',label = 'Pseudo-spectral',linewidth=1.5)
        if do_compare_with_FD:
            plt.plot(x,FD_V,color = 'blue',label = 'Finite-differences',marker='.',markevery=10)
        plt.ylim(-1,1)
        plt.legend()
        plt.show()

        # saves figure as pdf
        filename = 'figures/wave_plot_' + str(it) + '.pdf'
        plt.savefig(filename)
        print("plotted as ",filename)
        print("")
    
    # stability check
    if max(V) > 1.e3:
        print("simulation became unstable and blew up, please reduce time step size")
        sys.exit(1)

print("")
print("done")
print("")


## derivative of example functions to evaluate pseudo-spectral accuracy
if 1 == 1:
    print("derivative comparison example:")

    nx = 10000
    L = 2.0
    dx = L * np.pi * 2 / nx
    x = np.arange(0,L * np.pi * 2 * (nx-1)/nx + dx/2,dx)

    dk = 2.0 * np.pi / dx / nx

    print("  grid points nx = ",nx)
    print("  dx = ",dx," -> dk = ",dk)

    y = np.zeros(nx)

    # choose example: 1 == sin-function / 2 == step-function / 3 == triangle-function
    example_type = 1

    if example_type == 1:
        ## sin-function
        print("  example sin-function")

        y = np.sin(x);

    elif example_type == 2:
        ## step-function
        print("  example step-function")
        y[0:int(nx/3)] = 0; y[int(nx/3):2*int(nx/3)] = 1; y[2*int(nx/3):nx] = 0

    elif example_type == 3:
        ## triangle-function
        print("  example triangle-function")
        y[0:int(nx/2)] = np.arange(0,nx/2); y[int(nx/2):nx] = np.arange(int(nx/2),0,-1)

    # Fourier transform
    fft_y = np.fft.fft(y)

    # k indexing for numpy fft
    if nx % 2 == 0:
        k = np.asarray([i for i in range(0, int(nx/2))] + [0] + [i for i in range(- int(nx/2) + 1, 0)])
    else:
        k = np.asarray([i for i in range(0, int((nx-1)/2))] + [0] + [i for i in range(-int((nx-1)/2), 0)])

    # multiplication with i * k
    fft_y = 1j * k * dk * fft_y

    # inverse Fourier transform
    dy_FFT = np.real(np.fft.ifft(fft_y))

    # comparison with finite-difference solution (center difference)
    dy_FD = np.zeros(len(y))
    dy_FD[1:nx-1] = (y[2:nx]-y[0:nx-2]) / (2*dx)
    dy_FD[0] = (y[1]-y[0]) / dx
    dy_FD[nx-1] = (y[nx-1]-y[nx-2]) / dx

    ## sin-function example error
    if example_type == 1:
        # maximum errors
        ymax_FFT = max(dy_FFT-np.cos(x))
        ymax_FD = max(dy_FD-np.cos(x))
        print("  maximum error: FFT = ",ymax_FFT,"  FD = ",ymax_FD)

    # initial function plot
    plt.clf()
    plt.title('f(x)')
    plt.plot(x,y,color='black',label='f(x)')
    plt.legend()
    plt.show()

    # compares pseudo-spectral and finite-difference solution
    plt.clf()
    plt.title('derivative of f(x)')
    plt.plot(x,dy_FFT,color='red',label='FFT',marker='.',markevery=200)
    plt.plot(x,dy_FD, color='green',label='FD')
    # sine-function example
    if example_type == 1:
        plt.plot(x,np.cos(x),color='blue',label='analytical')

    plt.legend()
    plt.show()

    # sine-function example
    if example_type == 1:
        # errors in the approximation
        plt.clf()
        plt.title('errors')
        plt.plot(x,dy_FFT-np.cos(x),color='red',label='error FFT',marker='.',markevery=200)
        plt.plot(x,dy_FD-np.cos(x),color='green',label='error FD')
        ymax = max(ymax_FFT,ymax_FD)
        plt.ylim([-abs(ymax),abs(ymax)])
        plt.legend()
        plt.show()

