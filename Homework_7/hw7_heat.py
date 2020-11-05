#!/usr/bin/env python
#
# Homework 7: unsteady-state heat equation with SEM
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
uniform_elements = True
homogeneous_material = True

# model parameters
L     = 3000.0   # length (m)
rho_0 = 2500.0   # density (kg/m**3)
ca_0  = 300.0    # heat capacity ( cal/kg/K )
ka_1  = 1.0      # thermal conductivity ( cal/m/s/K )
ka_2  = 0.2

## meshing parameters
# number of spectral elements
NSPEC = 10

# number of GLL points (polynomial degree plus one)
NGLL = 3

################################################

# number of global points
NGLOB = (NGLL-1) * NSPEC + 1

# user output
print('SEM unsteady-state heat equation: ')
print('  number of spectral-elements: ',NSPEC)
print('  number of GLL points       : ',NGLL)
print('  number of global points    : ',NGLOB)
print('  Lagrange polynomial degree : ',NGLL - 1)
print('')

# estimate the time step 'time_step'
##>TODO: put your code here
DT = 5.e7
#total_dur = 76100.0  # 76.1 kyrs
total_dur = 10000.0    # 10.0 kyrs

# Julian astronomical year has 365.25 days
# 365.25 x 24 x 3600 s/hour = 31557600 seconds
num_julian_secs = 31557600.0

NSTEP = int(total_dur / (DT/num_julian_secs))

print("time step size: ",DT,"sec = ",str(np.around(DT/num_julian_secs,decimals = 2)),"yrs")
print("total duration: ",str(np.around(total_dur/1000.0,decimals = 2)),"kyrs")
print("number of time steps: ",NSTEP)

# start/end point location
X1 = 0; X2 = L

# regular grid
He = (X2-X1)/NSPEC
# point locations
x1 = np.zeros(NSPEC)
x2 = np.zeros(NSPEC)
if uniform_elements:
    # evenly spaced achors between 0 and 1
    print('even element spacing:')
    for ispec in range(NSPEC):
        x1[ispec] = X1 + ispec * He        # left point
        x2[ispec] = x1[ispec] + He         # right point
else:
    # unevenly spaced anchors between 0 and 1
    print('uneven element spacing:')
    left_length = L/3
    right_length = L - left_length
    nleft = int(NSPEC/2)
    nright = NSPEC - nleft
    for ispec in range(NSPEC):
        if ispec < nleft:
            x1[ispec] = ispec    * left_length/nleft
            x2[ispec] = (ispec+1)* left_length/nleft
        else:
            x1[ispec] = left_length + (ispec-nleft)  * right_length/nright
            x2[ispec] = left_length + (ispec-nleft+1)* right_length/nright
#debug
#for i in range(NSPEC): print('elem ',i+1,'size:',x2[i]-x1[i],'range:',x1[i],x2[i])

# element sizes
he = x2 - x1

# discretized material parameters
rho = rho_0 * np.ones([NGLL,NSPEC]) # density
ca = ca_0 * np.ones([NGLL,NSPEC])   # capacity
kappa = np.zeros([NGLL,NSPEC])      # conductivity
if homogeneous_material:
    # uniform conductivity
    for ispec in range(NSPEC):
        kappa[:,ispec] = ka_1 # constant conductivity
else:
    # non-uniform conductivity
    Nel2 = int(NSPEC/2)
    for ispec in range(NSPEC):
        if ispec < Nel2:
            kappa[:,ispec] = ka_1
        else:
            kappa[:,ispec] = ka_2

# Jacobian
dxidx = np.zeros([NGLL,NSPEC])
jacobian = np.zeros([NGLL,NSPEC])
for ispec in range(NSPEC):
    for i in range(NGLL):
        dxidx[i,ispec] = 2.0 / he[ispec]
        jacobian[i,ispec] = he[ispec] / 2.0

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
T_left  = ..    # temperature at point x = 0
T_right = ..     # temperature at point x = L

## initialization
K = np.zeros([NSPEC,NSPEC]) # stiffness matrix

d = np.zeros(NGLOB)     # temperature
d_dot = np.zeros(NGLOB) # d/dt temperature

## initial condition
d[:] = 0.0
# temperature on boundary
d[0] = T_left
d[NGLOB-1] = T_right

# figure plotting
fig,ax = plt.subplots(2,2,figsize = (8,8),sharex = True,sharey = True)
ax = ax.ravel()
isubplot = 0

# time loop
for itime in np.arange(0,NSTEP):
    # update temperature
    # predictor
    ##>TODO: put your code here
    d = ..


    ## solver
    rhs_global = np.zeros(NGLOB)
    # element loop
    for ispec in np.arange(0,NSPEC):

        # GLL point loop
        for i in range(NGLL):
            iglob = ibool[i,ispec]

            # local contribution
            ##>TODO: put your code here
            rhs_local = 0.0
            ..


            # boundary conditions
            ##>TODO: put your code here
            if fixed_BC:
                if (ispec == 0 and i == 0):
                    # left side
                    rhs_local = ..
                elif (ispec == NSPEC and i == NGLL-1):
                    # right side
                    rhs_local = ..

            # assembly
            ##>TODO: put your code here
            rhs_global[iglob] = ..

    # temperature increment
    # solve linear system: M d_dot + K d = F  for unknowns in d
    ##>TODO: put your code here
    d_dot[:] = ..

    # time scheme: corrector term
    ##>TODO: put your code here
    d[:] = ..

    # fixed temperature on boundary
    d[0] = T_left
    d[NGLOB-1] = T_right

    ## exact solution (for half-space)
    time = (itime+1)*DT
    # time in years
    time_yr = time / num_julian_secs
    N_ex = 200
    x_ex = np.linspace(X1,X2,N_ex)
    # exact solution for half-space problem
    T_ex = T_left * erfc(x_ex[:] / (2*np.sqrt(ka_1/rho_0/ca_0 * time)))

    # Plotting
    #dtime = int((NSTEP+1)/4.0) + 1
    dtime = int( np.linspace(1,NSTEP,4)[isubplot] ) - 1
    if itime == dtime:
        ax[isubplot].plot(x_ex,T_ex,label = 'Exact solution',color = 'blue')
        # figure x/y values
        figx = x.copy()
        figy = d.copy()
        ax[isubplot].plot(figx,figy,label = 'SEM solution',marker ='+',color = 'red',linewidth=0.85)
        # labels
        ax[isubplot].set_xlabel('x')
        ax[isubplot].set_ylabel('T')
        leg = ax[isubplot].legend(loc = 1)
        # titles
        ax[isubplot].set_title('time = '+ str(np.around(time_yr/1000.0,decimals = 1))+' kyrs')
        fig.suptitle('SEM solution - Homogeneous solution')
        isubplot += 1
        print("  done plot ",isubplot," out of 4")
        
plt.show()



