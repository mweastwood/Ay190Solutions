#!/usr/bin/env python

from __future__ import division
from pylab import *

f = lambda E,e,omega,t: E - omega*t - e*sin(E)
fprime = lambda E,e: 1 - e*cos(E)

def find_root(f,fprime,x0,N,max_step=1):
    '''
    This function uses Newton's method to find a root of f
    (with derivative fprime), given the initial guess x0,
    and the maximum number of iterations N.
    '''

    def newton_step(x):
        step = -f(x)/fprime(x)
        s = sign(step)
        step = s*minimum(abs(step),max_step)
        return x + step

    x = newton_step(x0)
    iterations = 1
    while iterations < N:
        if x > 1e-10 and (x-x0)/x < 1e-10: return (x,iterations)
        if x < 1e-10 and (x-x0)   < 1e-10: return (x,iterations)
        x0 = x
        x  = newton_step(x0)
        iterations += 1
    println("Warning: maximum number of iterations reached.")
    return (x,iterations)

# (a) Orbital parameters for the Earth
P = 365.25635
omega = 2*pi/P
e = 0.0167
a = 1 # AU
b = a*sqrt(1-e**2)

t = linspace(0,P)
roots = array([find_root(lambda x: f(x,e,omega,t_i),
                         lambda x: fprime(x,e), 0, 100) for t_i in t])
E = roots[:,0]; iterations = roots[:,1]
x1 = a*cos(E); y1 = b*sin(E)

print 'e = %.5f' % e
for t in [91,182,273]: # explicitly do the times called for in the question
    root = find_root(lambda x: f(x,e,omega,t), lambda x: fprime(x,e), 0, 100)
    print 't = %3d; E = %.3f; x = %+.3f; y = %+.3f; iterations = %d' \
                % (t,root[0],a*cos(root[0]),b*sin(root[0]),root[1])

# (b) Oops, the Earth got bumped a little bit
P = 365.25635
omega = 2*pi/P
e = 0.99999
a = 1 # AU
b = a*sqrt(1-e**2)

t = linspace(0,P)
roots = array([find_root(lambda x: f(x,e,omega,t_i),
                         lambda x: fprime(x,e), 0, 100) for t_i in t])
E = roots[:,0]; iterations = roots[:,1]
x2 = a*cos(E); y2 = b*sin(E)

print 'e = %.5f' % e
for t in [91,182,273]: # explicitly do the times called for in the question
    root = find_root(lambda x: f(x,e,omega,t), lambda x: fprime(x,e), 0, 100)
    print 't = %3d; E = %.3f; x = %+.3f; y = %+.3f; iterations = %d' \
                % (t,root[0],a*cos(root[0]),b*sin(root[0]),root[1])

# Plot the results
figure()
plot(x1,y1,'k-',label=r'$e=0.0167$')
plot(x2,y2,'r-',label=r'$e=0.99999$')
xlabel(r'$x/{\rm AU}$',fontsize=20)
ylabel(r'$y/{\rm AU}$',fontsize=20)
legend(loc='best')
gca().set_aspect('equal')
savefig('problem1.pdf')
close()

