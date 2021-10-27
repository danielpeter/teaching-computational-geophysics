#!/usr/bin/env python
#
# Homework 6: unsteady-state heat equation
#
# for conversion of matlab code to numpy, see:
# https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html
#
from __future__ import print_function
import os,sys

import numpy as np
print("")
print("numpy version: ",np.__version__)
print("")

#from scipy.linalg import solve
from scipy.special import erfc,erf

import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 20, 5
plt.rcParams['lines.linewidth'] = 0.5

## setup
# 1 - harmonic case,
# 2 - harmonic case, variable conductivity
# 3 - half-space cooling
problem_type = 1
print('problem: ',problem_type)

## meshing parameters
# number of elements
Nel = 10   # 10 or 20
# number of points
Np = Nel + 1

# start/end point location
if problem_type == 1 or problem_type == 2:
    X1 = 0; X2 = np.pi/2;
elif problem_type == 3:
    # half-space cooling
    X1 = 0; X2 = 20
else:
    print("problem type not recognized: ",problem_type)
    sys.exit(1)


# regular grid
He = (X2-X1)/Nel
# point locations
x1 = np.zeros(Nel)
x2 = np.zeros(Nel)
for i in range(Nel):
    x1[i] = X1 + i* He         # left point
    x2[i] = x1[i] + He         # right point

# irregular grid
if 1 == 0:
    He = 2.0/3.0 * (X2-X1)/Nel
    Nel2 = int(Nel/2)
    # point locations
    x1[0:Nel2] = X1 + np.arange(0,Nel2) * 2 * He
    x1[Nel2:Nel] = x1[Nel2-1] + 2 * He + np.arange(0,Nel2) * He
    x2[0:Nel2] = x1[0:Nel2] + 2 * He
    x2[Nel2:Nel] = x1[Nel2:Nel] + He
    #debug
    #for i in range(Nel): print('elem ',i+1,'size:',x2[i]-x1[i],'range:',x1[i],x2[i])

# element sizes
he = x2 - x1

## boundary condition
if problem_type == 1 or problem_type == 2:
    T1 = 1.0 # initial temperature at point x = L
    q0 = 0.0 # heat flux at point x = 0
elif problem_type == 3:
    # half-space cooling
    T0 = 0.0  # surface temp
    Tm = 1.0  # mantle temp
    q0 = 0.0  # heat flux at point x = 0

# discretized material parameters
rhoc = np.ones(Nel)
ka = np.ones(Nel)   #constant conductivity
if problem_type == 1:
    # uniform conductivity
    ka[0:Nel] = 1.0 # constant conductivity
elif problem_type == 2:
    # non-uniform conductivity
    Nel2 = int(Nel/2)
    ka[0:Nel2] = 2.0; ka[Nel2:Nel] = 1.0
elif problem_type == 3:
    # half-space cooling
    ka[0:Nel] = 1.0 # constant conductivity


# constant force factor
f_const = 0.0
f = f_const * np.ones(Np)

## time marching parameters
if problem_type == 3:
    # half-space cooling
    dt = 0.01 * min(he**2 * rhoc/ka)
    #dt = 0.016  # smaller dt leads to better solution
    Ntime = int(9.0/dt)
else:
    # harmonic case
    dt = 0.1 * min(he[0:Nel]**2 * rhoc[0:Nel]/ka[0:Nel])
    Ntime = 1000

alpha = 0.5

print("time step size: ",dt)
print("number of time steps: ",Ntime)
print("predictor-corrector scheme: ",alpha)

## initialization
M = np.zeros([Nel,Nel]) # mass (capacity) matrix
K = np.zeros([Nel,Nel]) # stiffness matrix
F = np.zeros(Nel)       # force (right-hand-side) vector

d = np.zeros(Nel)
d_dot = np.zeros(Nel)


## initial condition
if problem_type == 1 or problem_type == 2:
    # single harmonic initial
    d = d + np.cos(x1[0:Nel]) + 1
elif problem_type == 3:
    # half-space cooling
    d = Tm * np.ones(Nel)

# figure plotting
fig,ax = plt.subplots(2,2,figsize = (8,8),sharex = True,sharey = True)
ax = ax.ravel()
isubplot = 0

# time loop
for itime in np.arange(0,Ntime):
    ##> TODO: predictor
    d = ..
    d_dot = ..
    
    ## solver
    for e in np.arange(0,Nel):
        # number of local shape functions
        Nen = 2
        
        ## local to global view
        # setup local to global and equation numbering
        ID = np.arange(0,Nel+1)
        ID[Nel] = -1  # negative, index will get ignored in the assembly stage
        IEN = [e, e+1]

        ##> TODO: Location matrix: setup local to global numbering
        # (entries for this element)
        LM = np.empty(Nen,dtype=int)
        LM[:] = ..

        if itime == 0:
            print('Location matrix: ' + str(int(LM[0])) + '  ' + str(int(LM[1])))

        #> TODO: setup local mass, stiffness matrix and rhs vector
        me = ..
        ke = ..
        fe = ..

        # boundaries
        if e == 0:
            # left boundary
            fe[0] = ..
        elif e == Nel-1:
            # right boundary
            if problem_type == 1 or problem_type == 2:
                fe[0] = ..
            elif problem_type == 3:
                # half-space cooling
                fe[0] = ..
        
        ## assembly
        # add local contribution to the global matrix and rhs
        ind = np.nonzero(LM[:]+1) # note: nonzero(X) returns only the non-zero elements from vector X
        ind = ind[0]

        # check sub-matrix indexing np.ix_(..)
        # https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html
        #print(ind,LM[ind],M[np.ix_(LM[ind],LM[ind])])
        #for i in ind:
        #    iglob = LM[i]
        #    if itime == 0: print("e: ",e,"index i ",i, "iglob ",iglob,me[i,i],ke[i,i],fe[i])

        ##> TODO: construct global mass matrix
        M[np.ix_(LM[ind],LM[ind])] += ..
        # global stiffness matrix
        K[np.ix_(LM[ind],LM[ind])] += ..
        # global force vector
        F[np.ix_(LM[ind])] += ..

    # debug output
    #if itime == 0:
    #    print('time:',(itime+1)*dt)
    #    print(M)
    #    print(K)
    #    print(F)

    ##> TODO: solve linear system: M d_dot + K d = F  for unknowns in d
    R = ..
    d_dot = ..

    ##>TODO: corrector
    d = ..
    d_dot = ..


    ## exact solution
    # exact solution setup
    time = (itime+1)*dt
    N_ex = 200
    x_ex = np.linspace(X1,X2,N_ex)
    if problem_type == 1 or problem_type == 2:
        # harmonic
        T_ex = np.cos(x_ex[:]) * np.exp(-dt * itime) + T1
    elif problem_type == 3:
        # half-space cooling
        T_ex = erf((X2 - x_ex[:]) / (2*np.sqrt(ka[0]/rhoc[0] * time)))
    # Plotting
    dtime = int((Ntime+1)/4) + 1
    if itime%dtime == 0:
        # figure x/y values
        figx = x1.reshape(-1,1)
        figx = np.vstack([figx,X2]) # adds X2 to the end of array
        figy = d.reshape(-1,1)
        if problem_type == 1 or problem_type == 2:
            figy = np.vstack([figy,T1])  # adds T1 to the end of array
        elif problem_type == 3:
            figy = np.vstack([figy,T0])  # adds T0 to the end of array
        # plot
        ax[isubplot].plot(x_ex,T_ex,label = 'Exact solution',color = 'blue',linewidth=1.5)
        ax[isubplot].plot(figx,figy,label = 'FEM solution',marker ='o',color = 'red',linewidth=0.85)
        # labels
        ax[isubplot].set_xlabel('x')
        ax[isubplot].set_ylabel('T')
        leg = ax[isubplot].legend(loc = 3)
        # titles
        ax[isubplot].set_title('time = '+ str(np.around(time,decimals = 4)))
        if problem_type == 1:
            fig.suptitle('FEM solution - Homogeneous solution')
        elif problem_type == 2:
            fig.suptitle('FEM solution - Inhomogeneous solution')
        elif problem_type == 3:
            fig.suptitle('FEM solution - Half-space cooling solution')

        isubplot += 1
        
plt.show()

# saves figure as pdf
filename = 'figures/fem_solution.pdf'
plt.savefig(filename)
print("plotted as ",filename)
print("")

print("")
print("done")
print("")

