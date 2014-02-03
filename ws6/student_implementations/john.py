#!/usr/bin/env python

import numpy as np
from timeit import timeit
import matplotlib.pyplot as pl

# The square root of negative one is i, not j
i = 1j

# Function to compute a discrete Fourier Transform
def dft(x):
    N = len(x)
    w = np.exp(-(2*np.pi*i)/N)
    theMatrix = []
    for j in range(N):
        row = []
        for k in range (N):
            row.append(w ** (j*k))
        theMatrix.append(row)
    return np.dot(theMatrix, x)

# Abandoned attempt at faster version using map
#def dft2(x):
#    N = len(x)
#    w = np.exp(-(2*np.pi*i)/N)
#    theMatrix = []
#    for j in range(N):
#        k
#        row = map(w **)

