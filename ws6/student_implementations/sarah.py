#!/usr/bin/python/

import numpy as np
from timeit import timeit

def dft(x):
    # x,w
    N = x.shape[0]
    w = np.zeros((N,N),complex)

    # Fill w
    for k in range(0,N):
 	w[:,k] = np.exp(-np.linspace(0.,N-1,N)*2.*np.pi*1j*float(k)/float(N))

    # y
    y = np.dot(w,x)

    return y

