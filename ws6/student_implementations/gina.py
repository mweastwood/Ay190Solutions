import numpy as np
from scipy import pi
from timeit import timeit
import matplotlib.pyplot as plt

def dft(x):
    N = len(x)
    W=np.array([[np.exp(-2.0*pi*1j/N)]*N]*N)
    for k in range(N):
        for j in range(N):
            W[j,k]=(W[j,k])**(j*k)
    y=np.dot(W,x)
#    print x
#    print y
    return y

