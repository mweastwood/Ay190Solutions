#!/usr/bin/env python
import numpy as np
import scipy as sci

i=1j

def nyquist(h):
    return np.pi/h

def dft(func):
    N = len(func)
    w = np.exp(-2*np.pi*i/N)
    m = np.matrix([ [(w**j)**k for j in range(N)] for k in range(N)])
    func_hat = np.array(m.dot(func)).flatten()
    return func_hat

def fft_setup(func):
    N = len(func)
    power = int(np.ceil(np.log2(N)))
    if not power==np.log2(N):
        func = np.append(func,np.zeros(2**power-N))
        return func
    else:
        return func
                
def fft(func):
    func = fft_setup(func)
    N = len(func)
    func_hat = np.zeros(len(func)) + 0j
    for j in range(len(func)): # O(N)
        w = np.exp(-2*np.pi*i*j/N)
        func_hat[j] = fft_recursion(func,j,w)
    return func_hat
    
def fft_recursion(func, j, w):
    N = len(func)
    if N==1:
        return func[0]
    else: # Recursion should be O(log2(N)) 
        even = func[np.array([2*n     for n in range(N/2)])]
        odd  = func[np.array([2*n + 1 for n in range(N/2)])]
        return fft_recursion(even,j,w) + w*fft_recursion(odd,j,w)
