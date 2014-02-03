#!/usr/bin/env python
import numpy as np
from numpy import arange, array
from numpy import abs, exp, sin,cos,sqrt
import matplotlib.pyplot as pl

# Define the imaginary number i

i = 1j

# Problem a
def dft(x):
    N = len(x)
    w = exp(-2 * np.pi * i / N)
    A = []
    for j in range(N):
        row = []
        for k in range(N):
            row.append(w ** (j * k))
        A.append(array(row))
    return np.dot(array(A),array(x))

