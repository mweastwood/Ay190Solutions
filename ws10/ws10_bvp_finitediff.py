#!/usr/bin/env python

from __future__ import division
from pylab import *

# Boundary values
x_boundary = (0,1)
y_boundary = (0,0.1)

# ODE
f = lambda x: 12*x-4

# Setup the grid
N = 100 # Number of points that are not on the boundary
x = linspace(x_boundary[0],x_boundary[1],N+2)[1:-1]
h = x[1]-x[0]

# Construct the discrete Laplacian matrix
e = ones(N)
A = 1/h**2*(diag(e[1:],-1)-2*diag(e,0)+diag(e[1:],+1))

# Construct the RHS
b = f(x)
b[ 0] -= y_boundary[0]/h**2
b[-1] -= y_boundary[1]/h**2

# Solve
y = solve(A,b)

plot(x,y,'r-')
plot(x,2*x**3-2*x**2+0.1*x,'k-')
show()

