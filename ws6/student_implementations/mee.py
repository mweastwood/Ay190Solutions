# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 21:17:01 2014

Library for Ay190 Computational Astrophysics - WS 03

@author: Chatarin (Mee) Wong-u-railertkun
"""
#--------------------------------------------
import numpy as np

#--------------------------------------------
# Work Sheet 6
# Fourier Transform
#--------------------------------------------

def dft(x):
    n = len(x)
    
    a = np.asarray(range(n)*n).reshape(n,n)
    b = np.transpose(a)*1j
    
    c = np.mat(np.exp(a*b*-2.*np.pi/n))
    
    y = c*np.transpose(np.mat(x))
    return np.asarray(np.transpose(y)).reshape(n)