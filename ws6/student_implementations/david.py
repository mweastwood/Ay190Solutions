import pylab
from numpy import *
from matplotlib.pyplot import *
i=1.j
    
def dft(x):
    N = len(x)
    k=range(0,N)
    b=([(w(N)**k)**z for z in range(0,N)])
    return array([dot(b[s],x) for s in range(N)])

def w(N):
    return exp(-2*pi*i/N)

