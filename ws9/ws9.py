#!/usr/bin/env python

from __future__ import division
from pylab import *
import time

def gauss(A,b):
    A = copy(A); b = copy(b)
    N = A.shape[0]

    # Gaussian elimination
    for i in range(N):
        r = A[i,i]
        if abs(r) < 1e-13: print 'Warning: pivoting is not implemented'
        b[i]   /= r
        A[i,:] /= r
        for j in range(i+1,A.shape[0]):
            b[j]   -= A[j,i]*b[i]
            A[j,:] -= A[j,i]*A[i,:]

    # Back substitution
    x = zeros(N)
    for i in reversed(range(N)):
        if i == N:
            x[i] = b[i]
        else:
            x[i] = b[i]-sum(A[i,i+1:]*x[i+1:])
    return x

files = [('LSE1_m.dat','LSE1_bvec.dat'),
         ('LSE2_m.dat','LSE2_bvec.dat'),
         ('LSE3_m.dat','LSE3_bvec.dat'),
         ('LSE4_m.dat','LSE4_bvec.dat'),
         ('LSE5_m.dat','LSE5_bvec.dat')]

iterations = 10
for f in files:
    A = genfromtxt(f[0])
    b = genfromtxt(f[1])
    N,M = A.shape
    print 'Solving %d by %d system of linear equations...' % (N,M)

    t0 = time.clock()
    for i in range(iterations):
        gauss(A,b)
    t1 = time.clock()
    for i in range(iterations):
        linalg.solve(A,b)
    t2 = time.clock()

    print ' --> Gaussian elimination took %.4f seconds' % ((t1-t0)/iterations)
    print ' -->     LU decomposition took %.4f seconds' % ((t2-t1)/iterations)

