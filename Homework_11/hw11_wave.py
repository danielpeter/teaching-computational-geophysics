#!/usr/bin/env python
#
# Homework 11: wave equation with finite-volume method
#
#
from __future__ import print_function

import sys
import math
import numpy as np

from scipy.linalg import solve
from scipy.special import erfc,erf

import matplotlib.pyplot as plt
#plt.style.use('ggplot')
#plt.rcParams['figure.figsize'] = 20, 5
#plt.rcParams['lines.linewidth'] = 0.5

#############################################
## setup
# user parameters
apply_BC = True
BC_condition = 'Neumann'  # 'Dirichlet' or 'Neumann'

# model parameters
L    = 100.0     # length (m)
rho0 = 1.0       # density (kg/m**3)
mu0  = 1.0       # shear modulus

## meshing parameters
# number of (grid cell) elements
Nel = 100

# domain start/end point location
X_left = 0; X_right = L

# irregular mesh
is_irregular = False

# heterogeneous material
is_heterogeneous = False

################################################

# number of centroid points
Np = Nel

# number of local cell vertex points
N_local = Nel * 2

# user output
print('FVM wave equation: ')
print('  number of cell elements    : ',Nel)
print('  number of centroid points  : ',Np)
print('  number of local cell points: ',N_local)
print('')
print('  Boundary conditions        : ',apply_BC)
if apply_BC:
    print('  Boundary conditions        : ',BC_condition)
print('')
print('  irregular mesh             : ',is_irregular)
print('  heterogeneous material     : ',is_heterogeneous)
print('')

# estimate the time step 
DT = 0.01 
# duration
T_dur = 60.0
NSTEP = int(T_dur/DT)

print("time step size      : ",DT)
print("number of time steps: ",NSTEP)
print("total duration      : ",NSTEP*DT)
print("")


# meshing
x1 = np.zeros(Nel) # cell corner left point
x2 = np.zeros(Nel) # cell corner right point
if is_irregular:
    # irregular grid
    he = 2.0/3.0 * (X_right-X_left)/Nel
    Nel2 = int(Nel/2)
    # point locations
    # large elements
    x1[0:Nel2] = X_left + np.arange(0,Nel2) * 2 * he
    x2[0:Nel2] = x1[0:Nel2] + 2 * he
    # small elements
    x1[Nel2:Nel] = x1[Nel2-1] + 2 * he + np.arange(0,Nel2) * he
    x2[Nel2:Nel] = x1[Nel2:Nel] + he
else:
    # regular grid
    he = (X_right-X_left)/Nel
    # point locations
    for e in range(Nel):
        x1[e] = X_left + e * he        # left point
        x2[e] = x1[e] + he         # right point

# cell centroid / mid-point locations
xp = np.zeros(Nel)
for e in range(Nel):
    xp[e] = x1[e] + 0.5 * (x2[e] - x1[e])

# cell size (lengths)
Ve = np.zeros(Nel)
for e in range(Nel):
    Ve[e] = x2[e] - x1[e]

## material properties
# set up the local mesh properties (local level for each cell, left/right vertices)
rho = rho0 * np.ones([2,Nel]) # density
mu = mu0 * np.ones([2,Nel])   # shear modulus
# heterogeneous material
if is_heterogeneous:
    # assigns different material on last quarter of domain
    istart = 0
    while x1[istart] < 75.0: istart += 1
    print("material discontinuity at: ",x1[istart])
    print("")
    for e in np.arange(istart,Nel):
        ## wave speed c = sqrt(mu/rho) -> c = sqrt(2) * c0 = sqrt( 2 mu/rho) -> mu = 2 * mu0
        rho[:,e] = rho[:,e]
        mu[:,e] = 2 * mu[:,e]


## matrix initialization
M = np.zeros(Np)        # mass matrix
K = np.zeros((Np,Np))   # stiffness matrix
F = np.zeros(Np)        # right-hand-side vector

## initialization
displ = np.zeros(Np)     # displacement (defined at cell centroids)
veloc = np.zeros(Np)     # velocity
accel = np.zeros(Np)     # acceleration

# set up the boundary conditions
BC_displ_left = 0.0
BC_displ_right = 0.0
BC_grad_left = 0.0
BC_grad_right = 0.0

## initial condition
for i in range(Np):
    displ[i] = displ[i] + np.exp( -(xp[i]-50.0)**2 * 0.1 )

# calculate the mass matrix
##>TODO: put your code here
for e in range(Nel):
    # contribution in each cell: rho(x_p) * V = rho_average * V
    M[e] = ..

# calculate the stiffness matrix
##>TODO: put your code here
for e in np.arange(Nel):
    # contribution from cell element
    mu_A = mu[0,e]  # left cell vertex (at x_A)
    mu_Aplus1 = mu[1,e]  # right cell vertex (at x_A+1)

    # inverse 1/distance to neighbor centroids
    if e == 0:
        d_pminus1 = 0.0  # no neighbor on left, not defined
    else:
        d_pminus1 = 1.0/(xp[e]-xp[e-1])

    if e == Nel-1:
        d_pplus1 = 0.0  # no neighbor on right, not defined
    else:
        d_pplus1 = 1.0/(xp[e+1]-xp[e])

    ## stiffness term initialization
    ke_pminus1 = 0.0
    ke_p = 0.0
    ke_pplus1 = 0.0

    ## stiffness contributions, given centroid points P-1, P and P+1
    #>TODO: put your code here
    ke_pminus1 = .. 
    ke_p = ..
    ke_pplus1 = ..

    ## adds correction terms to have continuous stresses at cell faces between neighboring cells
    # shear moduli delta
    # discontinuity at left cell vertex (at x_A)
    if e == 0:
      mu_deltaA = 0.0 # no left neighbor
    else:
      mu_deltaA = mu[1,e-1] - mu_A
    # discontinuity at right cell vertex (at x_A+1)
    if e == Nel-1:
      mu_deltaAplus1 = 0.0 # no right neighbor
    else:
      mu_deltaAplus1 = mu[0,e+1] - mu_Aplus1
    # correction factors to take average stresses
    # (continuous across interface)
    ##>TODO: put your code here
    if np.abs(mu_deltaA) > 0.0:
        # to have 1/2(mu_A_left + mu_A_right) \partial_x s_A
        #         = [ mu_A_right + 1/2(mu_A_left - mu_A_right) ] \parital_x s_A
        #         = mu_A_right \parital_x s_A + 1/2 (mu_A_left - mu_A_right) \parital_x s_A
        #         = ke + 1/2 mu_deltaA \parital_x s_A
        ke_pminus1 += ..
        ke_p -= ..
    if np.abs(mu_deltaAplus1) > 0.0:
        # to have 1/2(mu_A+1_left + mu_A+1_right) \partial_x s_A+1
        #         = [ mu_A+1_left + 1/2(mu_A+1_right - mu_A+1_right) ] \parital_x s_A+1
        #         = mu_A+1_left \parital_x s_A+1 + 1/2 (mu_A+1_right - mu_A+1_right) \parital_x s_A+1
        #         = ke + 1/2 mu_deltaA+1 \parital_x s_A+1
        ke_p -= ..
        ke_pplus1 += ..

    ## boundary conditions
    if apply_BC:
        # first cell element
        if e == 0:
            # boundary on left side
            if BC_condition == 'Dirichlet':
                ## Dirichlet boundary
                dx1 = 1.0/(xp[0] - x1[0])
                ##>TODO: put your code here
                ke_pminus1 = 0.0                # not defined
                ke_p = ..
                ke_pplus1 = ..

            elif BC_condition == 'Neumann':
                ## Neumann boundary
                ##>TODO: put your code here
                ke_pminus1 = 0.0                # not defined
                ke_p = ..
                ke_pplus1 = ..

        # last cell element
        elif e == Nel-1:
            # boundary on right side
            if BC_condition == 'Dirichlet':
                ## Dirichlet boundary
                dx2 = 1.0/(x2[Nel-1] - xp[Nel-1])
                ##>TODO: put your code here
                ke_pminus1 = ..
                ke_p = ..
                ke_pplus1 = 0.0                 # not defined
            elif BC_condition == 'Neumann':
                ## Neumann boundary
                ##>TODO: put your code here
                ke_pminus1 = ..
                ke_p = ..
                ke_pplus1 = 0.0                 # not defined

    ## assembly
    # add contributions to the (global) matrix
    if e == 0:
        K[e,e] = ke_p
        K[e,e+1] = ke_pplus1
    elif e == Nel-1:
        K[e,e-1] = ke_pminus1
        K[e,e] = ke_p
    else:
        K[e,e-1] = ke_pminus1
        K[e,e] = ke_p
        K[e,e+1] = ke_pplus1

#debug
print("M :",M)
print("K :",K)
print("")

# figure plotting
NFIGS = 6
fig,ax = plt.subplots(NFIGS,1,figsize = (8,8),sharex = True,sharey = True)
ax = ax.ravel()
isubplot = 0

# time loop
for itime in np.arange(0,NSTEP):
    # Newmark time scheme: predictor
    # "predict" displacement, velocity, acceleration
    ##>TODO: put your code here
    displ[:] = ..
    veloc[:] = ..
    accel[:] = ..

    ## solver
    ## right-hand side calculation
    # F(i) = K(i,j) * d(j) matrix-vector multiplication
    F[:] = np.matmul(K[:,:],displ[:])

    # updates acceleration
    #>TODO: put your code here
    accel[:] = ..

    # Newmark time scheme: corrector
    # "correct" acceleration, velocity, displacement
    #>TODO: put your code here
    veloc[:] = ..

    #print("accel",accel)
    #print("veloc",veloc)
    #print("displ",displ)

    # Plotting
    dtime = int( np.linspace(1,NSTEP,NFIGS)[isubplot] ) - 1
    if itime == dtime:
        # figure x/y values
        data_x = xp.copy()
        data_y = displ.copy()
        #ax[isubplot].plot(data_x,data_y,label = 'FVM solution',marker ='o',fillstyle='none',markersize=1,
        #                  color = 'blue',linewidth=0.85,linestyle='dotted')
        ax[isubplot].plot(data_x,data_y,label = 'FVM solution',marker ='o',fillstyle='none',markersize=1,
                          color = 'blue',linewidth=0.85)

        # labels
        #ax[isubplot].set_xlabel('x')
        ax[isubplot].set_ylabel('displacement',fontsize=10)
        leg = ax[isubplot].legend(loc = 1)
        # titles
        ax[isubplot].set_title('time = '+ str(np.around((itime+1)*DT,decimals = 1)),fontsize=10)
        if is_heterogeneous:
            fig.suptitle('FVM solution - Heterogeneous solution')
        else:
            fig.suptitle('FVM solution - Homogeneous solution')
        isubplot += 1
        print("  done plot ",isubplot," out of ",NFIGS)

        # saves data file
        fname = 'figures/fvm_data_'+str(isubplot)+'.dat'
        np.savetxt(fname, np.array((data_x,data_y)).T,fmt='%10.8f')
        print("  written: ",fname)

# show result
plt.subplots_adjust(hspace=0.5)
plt.show()

# saves figure as pdf
filename = 'figures/fvm_solution.pdf'
fig.savefig(filename)
print("")
print("plotted as ",filename)
print("")

print("")
print("done")
print("")



