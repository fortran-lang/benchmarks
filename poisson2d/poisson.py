#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime


M = 100      
a = 1/(M-1)
target = 1e-8
analytical_tol = 5.0e-3

# Create arrays to hold potential values
phi = np.zeros([M,M],float)
phiprime = np.zeros([M,M],float)

#Analytical solution to the problem
sol = np.zeros([M,M],float)
begin = datetime.datetime.now()
def realsol(x,y):
    return y*(1-y)*x**3
for i in range(M):
    for j in range(M):
        sol[i,j] = realsol(a*i,a*j)

#Boundary conditions
for i in range(M):
    phi[M-1,i] = i*a*(1-i*a)
    phiprime[M-1,i] = i*a*(1-i*a)

#Source term
def ro(x,y):
    return -(6*x*y*(1-y) - 2*x**3)


# Main loop
delta = 1.0
while delta>target:

    # Calculate new values of the potential
    a2 = a**2
    for i in range(1,M-1):
        for j in range(1,M-1):

            phiprime[i,j] = (phi[i+1,j] + phi[i-1,j] \
                             + phi[i,j+1] + phi[i,j-1])/4 \
                    + a2*ro(i*a,j*a)/4																											

    # Calculate maximum difference from old values
    delta = np.max(abs(phi-phiprime))

    # Swap the two arrays around
    phi,phiprime = phiprime,phi
maxsol = np.max(abs(sol))
maxphi = np.max(abs(phi))
sol = sol/maxsol
phi = phi/maxphi
assert np.sum(np.abs(sol-phi))/(M*M) < analytical_tol
end = datetime.datetime.now()
dif = end-begin
print(dif.total_seconds())