#!/usr/bin/env python3


# -*- coding: utf-8 -*-
"""
Cythonized
"""
import numpy as np
cimport numpy as np
import matplotlib.pyplot as plt
cimport cython
from libc.math cimport pow
import time
import datetime


cdef double ro(double x,double y):
    cdef double res
    res = -(6.0*x*y*(1.0-y)-2.0*pow(x,3.0))
    return res
cdef double boundary(double y):
    return y*(1-y)
cdef double analytical_solution(double x, double y):
    return y*(1-y)*x**3
@cython.cdivision(True)
@cython.boundscheck(False)
cdef void run_calculation():
    cdef int M, counter
    cdef double V,target,a, analytical_tol
    M = 100         # Grid squares on a side
    target = 1e-8   # Target accuracy
    a = 1.0/(M-1)
    analytical_tol = 5.0e-3

    counter = 0
    cdef np.ndarray[np.double_t,ndim=2] phi = \
        np.zeros([M,M],dtype=np.double)

    cdef np.ndarray[np.double_t,ndim=2] phiprime = \
        np.zeros([M,M],np.double)
    cdef np.ndarray[np.double_t,ndim=2] roa = \
        np.zeros([M,M],np.double)
    cdef np.ndarray[np.double_t,ndim=2] analytical_sol = \
        np.zeros([M,M],np.double)
    cdef int i,j
    cdef double delta
    cdef double a2
    a2 = a**2
    delta = 1.0

    for i in range(M):
        for j in range(M):
            roa[i,j] = ro(i*a, j*a)
    for i in range(M):
        phiprime[M-1,i] = boundary(i*a)
        phi[M-1,i] = boundary(i*a)
    for i in range(M):
        for j in range(M):
            analytical_sol[i,j] = analytical_solution(i*a,j*a)

    while delta>target:
        counter += 1
        for i in range(1,M-1):
            for j in range(1,M-1):
                phiprime[i,j] = (phi[i+1,j] + phi[i-1,j] \
                                 + phi[i,j+1] + phi[i,j-1])/4 \
                        + a2*roa[i,j]/4.0

    # Calculate maximum difference from old values
        delta = np.max(abs(phi-phiprime))
    # Swap the two arrays around
        phi,phiprime = phiprime,phi
    assert np.sum(np.abs(analytical_sol-phi))/(M*M) < analytical_tol
    
    print(counter)
    
begin = datetime.datetime.now()
run_calculation()
end = datetime.datetime.now()
dif = end-begin
print(dif.total_seconds())
