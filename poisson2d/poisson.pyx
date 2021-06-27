#!/usr/bin/env python3


# -*- coding: utf-8 -*-
"""
Cythonized
"""
import numpy as np
cimport numpy as np
import matplotlib.pyplot as plt
cimport cython
import time
import datetime


cdef int ro(double x,double y):
    cdef double startone, endone, starttwo, endtwo
    startone = 0.6
    endone = 0.8
    starttwo = 0.2
    endtwo = 0.4
    if x>startone and x<endone and y>startone and y<endone:
        return 1
    elif x>starttwo and x<endtwo and y>starttwo and y<endtwo:
        return -1
    else:
        return 0


@cython.cdivision(True)
@cython.boundscheck(False)
cdef void run_calculation():
    cdef int M, counter
    cdef double V,target,a,epsilon0
    M = 300         # Grid squares on a side
    V = 1.0         # Voltage at top wall
    target = 1e-6   # Target accuracy
    a = 0.01
    epsilon0 = 8.85e-12

    counter = 0
    cdef np.ndarray[np.double_t,ndim=2] phi = \
        np.zeros([M+1,M+1],dtype=np.double)

    cdef np.ndarray[np.double_t,ndim=2] phiprime = \
        np.zeros([M+1,M+1],np.double)
    cdef int i,j
    cdef double delta
    cdef double a2
    a2 = a**2
    delta = 1.0
    while delta>target:
        counter += 1
        for i in range(1,M):
            for j in range(1,M):
                phiprime[i,j] = (phi[i+1,j] + phi[i-1,j] \
                                 + phi[i,j+1] + phi[i,j-1])/4 \
                        + a2/(4*epsilon0)*ro(i*a,j*a)

    # Calculate maximum difference from old values
        delta = np.max(np.abs(phi-phiprime))

    # Swap the two arrays around
        phi,phiprime = phiprime,phi
    print(counter)
    
begin = datetime.datetime.now()
run_calculation()
end = datetime.datetime.now()
dif = end-begin
print(dif.total_seconds())
