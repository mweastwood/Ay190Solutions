#!/usr/bin/env python
import sys,math
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp


def apply_bcs(x,y):
    y[0] = y[1]
    y[len(y)-1] = y[len(y)-2]
    return y

def rhs(q,dx):
    n = len(q)
    dqdt = np.zeros(n)
    for i in range(1,n-1):
        if(q[i] > 0):
            dqdt[i] = -q[i]/(dx) * (q[i] - q[i-1])
        else:
            dqdt[i] = -q[i]/(dx) * (q[i+1] - q[i])
    return dqdt

def analytic(x,t):
    return 0.125*np.sin( 2*np.pi*x/100. ) 
    #return 0.0625+0.000625*(100-(x-30.)**2) * (x > 20) * (x < 40)


x = np.linspace(0,100,101)
dx = x[1]-x[0]

n = len(x)
y = np.zeros(n)
dt = 0.25 * dx / 0.125


#time integration loop
ntmax = 200
t = 0.0

#set up initial conditions
#for i in range(n):
y = analytic(x,t)


# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x,y,'x-')
mpl.show()

yold2 = y
yold = y
for it in range(ntmax):
    yold2 = yold
    yold = y
    k1 = dt*rhs(yold,dx)
    k2 = dt*rhs(yold+0.5*k1,dx)
    y = yold + k2
    y = apply_bcs(x,y)
    t = t + dt
    print "it = ",it," t = ",it*dt
    if it % 10 == 0:
        mpl.clf()
        mpl.plot(x,y,'x-')
        mpl.draw()

mpl.ioff()
mpl.clf()
mpl.plot(x,y)
mpl.show()


