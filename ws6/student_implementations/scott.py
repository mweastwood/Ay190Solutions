import numpy as np
from timeit import timeit
import matplotlib.pyplot as plt

def dft(x):
    N=len(x)
    W=np.zeros((N,N),dtype=complex)
#Calculate the W matrix
    for k in range(N):
        for j in range(N):
            W[k][j]=np.exp((-2.*np.pi*1J*j*k)/N)
#and dot it with x
    y=np.dot(W,x)
    return y

