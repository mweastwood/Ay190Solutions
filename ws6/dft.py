#!/usr/bin/env python

from __future__ import division
from pylab import *

def dft(x):
    N = len(x)
    w = exp(-2j*pi/N)
    A = zeros((N,N),dtype=complex)
    for i in range(N):
        for j in range(N):
            A[i,j] = w**(i*j)
    return dot(A,x)

