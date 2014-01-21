#!/usr/bin/env python

from __future__ import division
from pylab import *

# Define f and its derivative
f      = lambda x:   x**3 -  5*x**2 + x
fprime = lambda x: 3*x**2 - 10*x    + 1

def forward_difference(x,y):
    xprime = x[:-1]                        # Grid points where the derivative is estimated
    yprime = (y[1:]-y[:-1])/(x[1:]-x[:-1]) # Estimate of the derivative
    return (xprime,yprime)

def central_difference(x,y):
    xprime = x[1:-1]                       # Grid points where the derivative is estimated
    yprime = (y[2:]-y[:-2])/(x[2:]-x[:-2]) # Estimate of the derivative
    return (xprime,yprime)

def estimate_convergence_order(function,derivative,estimate):
    '''
    This function estimates the convergence order of the
    estimate of the derivative.
    - function: the function whose derivative we are evaluating
    - derivative: the analytic derivative of the function
    - estimate: a function that estimates the derivative
    '''
    grid1 = linspace(-2,6,101)
    grid2 = linspace(-2,6,201)

    # Estimate the derivative on both grids
    x1,y1 = estimate(grid1,function(grid1))
    x2,y2 = estimate(grid2,function(grid2))
    
    # We want to compare the derivatives at the same values of x,
    # so throw away all the values evaluated at different values
    y2 = array([y2[i] for i in range(len(y2)) if any(abs(x2[i]-x1)<1e-10)])
    x2 = array([x2[i] for i in range(len(x2)) if any(abs(x2[i]-x1)<1e-10)])

    convergence_factor = mean(abs((y1-derivative(x1))/(y2-derivative(x2))))
    convergence_order  = log(convergence_factor)/log(2)
    return convergence_order

print 'Forward differencing is n=%.1f convergent' \
        % estimate_convergence_order(f,fprime,forward_difference)
print 'Central differencing is n=%.1f convergent' \
        % estimate_convergence_order(f,fprime,central_difference)

