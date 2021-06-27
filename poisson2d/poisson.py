#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime

# Constants
M = 100         # Grid squares on a side
V = 1.0         # Voltage at top wall
target = 1e-6   # Target accuracy
a = 0.01
epsilon0 = 8.85e-12

# Create arrays to hold potential values
phi = np.zeros([M+1,M+1],float)
#phi[0,:] = V
phiprime = np.zeros([M+1,M+1],float)

def ro(x,y):
    if x>0.6 and x<0.8 and y>0.6 and y<0.8:
        return 1
    elif x>0.2 and x<0.4 and y>0.2 and y<0.4:
        return -1
    else:
        return 0


# Main loop
begin = datetime.datetime.now()
delta = 1.0
while delta>target:

    # Calculate new values of the potential
    a2 = a**2
    for i in range(1,M):
        for j in range(1,M):

            phiprime[i,j] = (phi[i+1,j] + phi[i-1,j] \
                             + phi[i,j+1] + phi[i,j-1])/4 \
                    + a2/4/epsilon0*ro(i*a,j*a)																												

    # Calculate maximum difference from old values
    delta = np.max(abs(phi-phiprime))

    # Swap the two arrays around
    phi,phiprime = phiprime,phi
end = datetime.datetime.now()
dif = end-begin
print(dif.total_seconds())
# Make a plot
plt.imshow(phi,origin='lower')
plt.savefig("purepython.png")
