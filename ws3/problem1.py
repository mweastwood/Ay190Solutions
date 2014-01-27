#!/usr/bin/env python

from __future__ import division
from pylab import *

def midpoint_rule(function,a,b,h):
    x = linspace(a,b,(b-a)/h+1)
    midpoint = function((x[1:]+x[:-1])/2)
    Q = (x[1:]-x[:-1])*midpoint
    return sum(Q)

def trapezoidal_rule(function,a,b,h):
    x = linspace(a,b,(b-a)/h+1)
    y = function(x)
    Q = 0.5*(x[1:]-x[:-1])*(y[1:]+y[:-1])
    return sum(Q)

def simpsons_rule(function,a,b,h):
    x = linspace(a,b,(b-a)/h+1)
    y = function(x)
    midpoint = function((x[1:]+x[:-1])/2)
    Q = (x[1:]-x[:-1])/6*(y[1:]+4*midpoint+y[:-1])
    return sum(Q)

N = floor(logspace(1,3,1e2)) # number of points
h = pi/(N-1) # step size

# (a) Integrate sin(x) using the trapezoidal rule and Simpson's rule
midp = array([   midpoint_rule(sin,0,pi,h_i) for h_i in h])
trap = array([trapezoidal_rule(sin,0,pi,h_i) for h_i in h])
simp = array([   simpsons_rule(sin,0,pi,h_i) for h_i in h])

figure()
plot(log10(h),log10(abs(midp-2)),'k-',label='midpoint rule')
plot(log10(h),log10(abs(trap-2)),'r-',label='trapezoidal rule')
plot(log10(h),log10(abs(simp-2)),'b-',label='Simpson\'s rule')
xlabel(r'$\log_{10}\,h$',fontsize=20)
ylabel(r'$\log_{10}\,\left|\int_0^\pi\sin\,x\,{\rm d}x - 2\right|$',fontsize=20)
xlim(log10(h.min()),log10(h.max()))
gca().invert_xaxis()
legend(loc='best')
savefig('problem1a.pdf')
close()

# (b) Integrate x*sin(x) using the trapezoidal rule and Simpson's rule
midp = array([   midpoint_rule(lambda x: x*sin(x),0,pi,h_i) for h_i in h])
trap = array([trapezoidal_rule(lambda x: x*sin(x),0,pi,h_i) for h_i in h])
simp = array([   simpsons_rule(lambda x: x*sin(x),0,pi,h_i) for h_i in h])

plot(log10(h),log10(abs(midp-pi)),'k-',label='midpoint rule')
plot(log10(h),log10(abs(trap-pi)),'r-',label='trapezoidal rule')
plot(log10(h),log10(abs(simp-pi)),'b-',label='Simpson\'s rule')
xlabel(r'$\log_{10}\,h$',fontsize=20)
ylabel(r'$\log_{10}\,\left|\int_0^\pi x\sin\,x\,{\rm d}x - \pi\right|$',fontsize=20)
xlim(log10(h.min()),log10(h.max()))
gca().invert_xaxis()
legend(loc='best')
savefig('problem1b.pdf')
close()

