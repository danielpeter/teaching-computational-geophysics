#!/usr/bin/env python
#
# Homework 8: wave equation with SEM
#
#
from __future__ import print_function

import sys
import math
import numpy as np

from scipy.linalg import solve
from scipy.special import erfc,erf

import matplotlib.pyplot as plt

#############################################
## setup
# user parameters
fixed_BC = True

# model parameters
L    = 100.0   # length (m)
rho0 = 1.0       # density (kg/m**3)
mu0  = 1.0       # shear modulus

## meshing parameters
# number of spectral elements
NSPEC = 20

# number of GLL points (polynomial degree plus one)
NGLL = 7

################################################

# number of global points
NGLOB = (NGLL-1) * NSPEC + 1

# user output
print('SEM wave equation: ')
print('  number of spectral-elements: ',NSPEC)
print('  number of GLL points       : ',NGLL)
print('  number of global points    : ',NGLOB)
print('  Lagrange polynomial degree : ',NGLL - 1)
print('')

# estimate the time step 'time_step'
##>TODO: put your code here
DT = 0.25
NSTEP = 220

print("time step size      : ",DT)
print("number of time steps: ",NSTEP)
print("total duration      : ",NSTEP*DT)
print("")

# start/end point location
X1 = 0; X2 = L
# evenly spaced anchors between 0 and 1
He = (X2-X1)/NSPEC
# point locations
x1 = np.zeros(NSPEC)
x2 = np.zeros(NSPEC)
print('even element size:',He)
for ispec in range(NSPEC):
    x1[ispec] = X1 + ispec * He        # left point
    x2[ispec] = x1[ispec] + He         # right point
# element sizes
he = x2 - x1

# set up the mesh properties
rho = rho0 * np.ones([NGLL,NSPEC]) # density
mu = mu0 * np.ones([NGLL,NSPEC])   # shear modulus

# Jacobian
dxidx = np.zeros([NGLL,NSPEC])
jacobian = np.zeros([NGLL,NSPEC])
for ispec in range(NSPEC):
    for i in range(NGLL):
        dxidx[i,ispec] = 2.0/he[ispec]
        jacobian[i,ispec] = he[ispec]/2.0

# set up local to global numbering
ibool = np.empty([NGLL,NSPEC],dtype=int)
iglob = 0
for ispec in range(NSPEC):
    for i in range(NGLL):
        if i >= 1: iglob = iglob + 1
        ibool[i,ispec] = iglob
#debug
for ispec in range (NSPEC): print('Location matrix: ',ibool[:,ispec])


# set up of the Gauss-Lobatto-Legendre points:
#   xigll  - GLL nodal positions in reference domain
#   wgll   - associated nodal weights for quadrature rule
#   hprime - derivatives of the Lagrange polynomials
#            hprime(i,j) = h'_i(xigll_j) by definition of the derivative matrix
#
# (values taken from output of fortran xdiffusion code)
if NGLL == 3:
    # second-order elements
    xigll = np.array([-1.0, 0.0, 1.0])
    wgll  = np.array([0.33333333333333333,1.3333333333333333,0.33333333333333333])
    hprime = np.array([ [-1.5,-0.5, 0.5], [2.0, 0.0,-2.0], [-0.5, 0.5, 1.5] ])
elif NGLL == 4:
    # third-order elements
    xigll = np.array([-1.0, -0.44721359549995798, 0.44721359549995798, 1.0])
    wgll  = np.array([0.16666666666666666,0.83333333333333333,0.83333333333333333,0.16666666666666666])
    hprime = np.array([ [-3.0, -0.80901699437494745, 0.30901699437494740, -0.5],\
                        [4.0450849718747381, 0.0, -1.1180339887498947, 1.5450849718747375],\
                        [-1.5450849718747375, 1.1180339887498947, 0.0, -4.0450849718747381],\
                        [0.5, -0.30901699437494740, 0.80901699437494745, 3.0] ])
elif NGLL == 5:
    # 4th-order elements
    xigll = np.array([-1.0, -0.65465367070797720, 0.0, 0.65465367070797709, 1.0])
    wgll  = np.array([0.1, 0.5444444444444444, 0.711111111111111, 0.5444444444444444, 0.1])
    hprime = np.array([ [-5.0,               -1.2409902530309833,   0.375,              -0.25900974696901718,  0.5],\
                        [6.7565024887242409,  0.0,                 -1.3365845776954530,  0.76376261582597338, -1.4101641779424265],\
                        [-2.6666666666666666, 1.7457431218879389,   0.0,                 -1.7457431218879396,  2.6666666666666666],\
                        [1.4101641779424268, -0.76376261582597338,  1.3365845776954532,  0.0,                 -6.7565024887242382],\
                        [-0.5,                0.25900974696901713, -0.375,               1.2409902530309824,   5.0000000000000000] ])
elif NGLL == 7:
    # 6th-order elements
    xigll = np.array([-1.0, -0.830223896278566964, -0.468848793470714231, 0.0, 0.468848793470714231, 0.830223896278566964, 1.0])
    wgll  = np.array([ \
      0.047619047619047616,  0.276826047361566130,  0.431745381209862611,0.487619047619047619,  0.431745381209862611,  0.276826047361566130,  0.047619047619047616 \
                     ])
    hprime = np.array([ \
      [-10.500000000000000000, -2.442926014244291455,  0.625256665515342092, -0.312500000000000000,  0.226099400942574691, -0.226611870395445392,  0.500000000000000000],\
      [14.201576602919812942,  0.000000000000000000, -2.215804283169970468,  0.907544471268820652, -0.616390835517579339,  0.602247179635785668, -1.317373435702434037],\
      [-5.668985225545506879,  3.455828214294285328,  0.000000000000000000, -2.006969240588753145,  1.066441904006374619, -0.961339797288711884,  2.049964813076742498],\
      [ 3.200000000000000178, -1.598606688098367368,  2.266698087085999180,  0.000000000000000000, -2.266698087085999180,  1.598606688098367368, -3.200000000000000178],\
      [-2.049964813076742498,  0.961339797288711884, -1.066441904006374619,  2.006969240588753145,  0.000000000000000000, -3.455828214294285328,  5.668985225545506879],\
      [ 1.317373435702434037, -0.602247179635785668,  0.616390835517579339, -0.907544471268820652,  2.215804283169970468,  0.000000000000000000,-14.201576602919812942],\
      [-0.500000000000000000,  0.226611870395445392, -0.226099400942574691,  0.312500000000000000, -0.625256665515342092,  2.442926014244291455, 10.500000000000000000]])
else:
    print("Invalid setting for NGLL ",NGLL," (only supported for NGLL == 3, 4, 5 or 7)")
    sys.exit(1)

# get the global grid points
x = np.empty(NGLOB)
for ispec in range(NSPEC):
    for i in range(NGLL):
        iglob = ibool[i,ispec]
        x[iglob] = 0.5 * (1.0-xigll[i])*x1[ispec] + 0.5 * (1.0+xigll[i])*x2[ispec]


# calculate the global mass matrix 'mass_global'
##>TODO: put your code here
M_global = np.zeros(NGLOB)
..

# set up the boundary conditions
##>TODO: put your code here
displ_1 = ..
displ_NGLOB = ..
grad_1 = ..
grad_NGLOB = ..

## initialization
displ = np.zeros(NGLOB)     # displacement
veloc = np.zeros(NGLOB)     # velocity
accel = np.zeros(NGLOB)     # acceleration

## add source
## initial condition
#>TODO: put your code here
..


# figure plotting
fig,ax = plt.subplots(4,1,figsize = (8,8),sharex = True,sharey = True)
ax = ax.ravel()
isubplot = 0

# time loop
for itime in np.arange(0,NSTEP):

    # "predict" displacement, velocity, initialize acceleration
    # Newmark time scheme: predictor terms
    # predictor
    ##>TODO: put your code here
    displ[:] = ..
    veloc[:] = ..
    accel[:] = ..

    # boundary conditions: Dirichlet boundary condition
    # (zero displacement at boundary points)
    if fixed_BC:
        displ[0] = displ_1
        displ[NGLOB-1] = displ_NGLOB

    ## solver
    # element loop
    for ispec in np.arange(0,NSPEC):

        # ! contribution from local element & assembly on global points
        #>TODO: put your code here
        ..

        # boundary conditions
        ##>TODO: put your code here
        ..

    # updates acceleration
    #>TODO: put your code here
    accel[:] = ..

    # "correct" acceleration, velocity, displacement
    # Newmark time scheme: corrector term
    #>TODO: put your code here
    veloc[:] = ..


    # Plotting
    dtime = int( np.linspace(1,NSTEP,4)[isubplot] ) - 1
    if itime == dtime:
        # figure x/y values
        figx = x.copy()
        figy = displ.copy()
        ax[isubplot].plot(figx,figy,label = 'SEM solution',marker ='+',color = 'blue',linewidth=0.85)
        # labels
        #ax[isubplot].set_xlabel('x')
        ax[isubplot].set_ylabel('displacement')
        leg = ax[isubplot].legend(loc = 1)
        # titles
        ax[isubplot].set_title('time = '+ str(np.around((itime+1)*DT,decimals = 1)))
        fig.suptitle('SEM solution - Homogeneous solution')
        isubplot += 1
        print("  done plot ",isubplot," out of 4")
        
plt.show()



