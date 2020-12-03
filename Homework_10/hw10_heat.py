#!/usr/bin/env python
#
# Homework 10: steady-state heat equation with finite-volume method
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

#############################################
## setup

## meshing parameters
# number of elements
Nel = 10     # 10 or 20

# domain start/end point location
X_left = 0.0; X_right = 1.0

#############################################

# meshing
x1 = np.zeros(Nel) # cell corner left point
x2 = np.zeros(Nel) # cell corner right point
if True:
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
        x1[e] = X_left + i * he        # left point
        x2[e] = x1[e] + he         # right point

# cell centroid / mid-point locations
xp = np.zeros(Nel)
for e in range(Nel):
    xp[e] = x1[e] + 0.5 * (x2[e] - x1[e])

# cell size (lengths)
Ve = np.zeros(Nel)
for e in range(Nel):
    Ve[e] = x2[e] - x1[e]


#debug
#for e in range(Nel): print('elem ',e+1,'size:',x2[e]-x1[e],'range:',x1[e],x2[e],'\tcentroid:',xp[e])

#> TODO: change boundary conditions
## boundary condition
T1 = 1            # Dirichlet: fixed temperature at point x = 1
q0 = 1            # Neumann: heat flux at point x = 0

# constant force factor
f_const = 5       # 0 or 5

# number of centroid points
Np = Nel

# given force
f = f_const * np.ones(Np);

## initialization
K = np.zeros((Np,Np))   # stiffness matrix
F = np.zeros(Np)        # rhs/force vector

print("")
print("number of elements = ",Nel)
print("")

# sets up linear system
for e in np.arange(0,Nel):

    #> TODO: element contributions and rhs vector
    ## rhs/force contribution
    fe = ..

    ## adds boundary conditions to rhs
    # left/right boundaries
    if e == 0:
        # left boundary: Neumann
        fe = ..
    elif e == Nel-1:
        # right boundary: Dirichlet
        fe = ..

    # (global) rhs vector
    F[e] = fe

    ## diffusion term
    ke_pminus1 = 0.0
    ke_p = 0.0
    ke_pplus1 = 0.0
    if e == 0:
        # Neumann boundary on left side
        ke_pminus1 = ..
        ke_p = ..
        ke_pplus1 = ..
    elif e == Nel-1:
        # Dirichlet boundary on right side
        ke_pminus1 = ..
        ke_p = ..
        ke_pplus1 = ..
    else:
        # neighbor contributions
        ke_pminus1 = ..
        ke_p = ..
        ke_pplus1 = ..

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

print("")
print("global stiffness K: \n",K)
print("global rhs F      : \n",F)
print("")
#print("inverted stiffness K^-1: \n",np.linalg.inv(K),"\n")

#> TODO: solve linear system: K d = F  for unknowns in d
d = ..

## exact solution
N = 200
x_ex = np.linspace(X_left,X_right,N)
T_ex = T1 + (1-x_ex)*q0 + (1-x_ex**2)*f_const/2

## plot result
plt.clf()
plt.title('Nel = ' + str(Nel) + ' f = ' + str(f_const))
# temperature values at cell centroid locations
plt.plot(xp,d,color = 'red',label='FVM',linewidth=1.5,marker='.')
# exact solution
plt.plot(x_ex,T_ex,color = 'blue',label='exact solution')
plt.legend()
fig = plt.gcf()

# show
plt.show()

# saves figure as pdf
filename = 'figures/fvm_solution.pdf'
fig.savefig(filename)

print("plotted as ",filename)
print("")

print("")
print("done")
print("")
