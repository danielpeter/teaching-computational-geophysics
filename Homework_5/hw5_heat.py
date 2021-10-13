#!/usr/bin/env python
#
# Homework 5: steady-state heat equation
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

## meshing parameters
# number of elements
Nel = 10     # 2 or 10

# start/end point location
x1 = 0.0; x2 = 1.0

# element size
he = (x2-x1)/Nel

# number of points
Np = Nel + 1
# point locations
x = np.linspace(x1,x2,Np)

## boundary condition
T1 = 1           # initial temperature at point x = 1
q0 = 1           # heat flux at point x = 0

# constant force factor
f_const = 1.0    # 0 or 1

# given force
f = f_const * np.ones(Np);

## initialization
K = np.zeros((Nel,Nel))   # stiffness matrix
F = np.zeros(Nel)         # force vector

print("")
print("number of elements = ",Nel)
print("")

# sets up linear system
for e in np.arange(0,Nel):
    ## number of local shape functions
    Nen = 2

    ## local to global view
    # sets up global to equation numbering
    ID = np.zeros(Np,dtype=int)
    ID[0:Nel] = np.arange(1,Nel+1)
    ID[Np-1] = 0   # no entry

    IEN = np.zeros(Nen,dtype=int)
    IEN[0:Nen] = [e,e+1]

    #> TODO: Location matrix: setup local to global numbering
    # (entries for this element)
    LM = np.zeros(Nen,dtype=int)
    LM[0:Nen] = ..
    print("Location matrix: ",LM)

    #> TODO: setup local stiffness matrix and rhs vector
    ke = ..
    fe = ..

    ## assembly
    # (add local contribution to the global matrix)
    ind = LM.nonzero()  # note: nonzero() returns only the non-zero elements from vector X
    indices = ind[0]    # for example indices can be : [0 1] or just [0]
    #debug
    #print("indices = ",indices," len = ",len(indices))
    # local index range; gets first and last entry from indices
    k1 = indices[0]; k2 = indices[-1]
    # global index range; -1 to have index range between 0 to n-1 for python arrays
    i1 = LM[k1]-1; i2 = LM[k2]-1

    #> TODO: construct global stiffness matrix
    K[i1:i2+1,i1:i2+1] = ..

    # global force vector
    F[LM[ind]-1] = ..

print("")
print("global stiffness K: \n",K)
print("global rhs F      : \n",F)
print("")

#> TODO: solve linear system: K d = F  for unknowns in d
d = ..

## exact solution
N = 200
x_ex = np.linspace(x1,x2,N)
T_ex = T1 + (1-x_ex)*q0 + (1-x_ex**2)*f_const/2

## plot result
plt.clf()
plt.title('Nel = ' + str(Nel) + ' f = ' + str(f_const))
plt.plot(x,np.concatenate((d,[T1])),color = 'red',label='FEM',linewidth=1.5,marker='.')
plt.plot(x_ex,T_ex,color = 'blue',label='exact solution')
plt.legend()
plt.show()

# saves figure as pdf
filename = 'figures/fem_solution.pdf'
plt.savefig(filename)
print("plotted as ",filename)
print("")

print("")
print("done")
print("")
