#!/usr/bin/env python

import numpy as np

def dft(x):

    N = len(x)
    w = np.exp(-2.*np.pi*1.j/N)

    grid = np.indices((N, N))
    ind = grid[0]*grid[1]

    wjj = np.power(w,ind)

    return np.dot(wjj, x)

