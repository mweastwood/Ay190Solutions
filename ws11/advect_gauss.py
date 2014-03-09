import sys,math
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp

def apply_bcs(x,y):
    # apply boundary conditions
    # you need to fill in code
    y[ 0] = y[ 1]
    y[-1] = y[-2]
    return y

def analytic(x,x0,sigma,v,t):
    return np.exp(-(x-x0-v*t)**2/(2.*sigma**2))

def rhs(x,y,v,derivative="central"):
    dydt = np.zeros(y.shape)
    if   derivative == "right":
        dydx = (y[1:]-y[:-1])/(x[1:]-x[:-1])
        dydt[:-1] = -v*dydx
    elif derivative == "left":
        dydx = (y[1:]-y[:-1])/(x[1:]-x[:-1])
        dydt[1:] = -v*dydx
    elif derivative == "central":
        dydx = (y[2:]-y[:-2])/(x[2:]-x[:-2])
        dydt[1:-1] = -v*dydx
    else:
        error("Unknown derivative type.")
    return dydt

# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x = np.linspace(0,100,1001)
# parameters
dx = x[1]-x[0]
v = 0.1

n = len(x)
y = np.zeros(n)
cfl = 1.0
dt = cfl*v/dx
t = 0.0

# for initial data
sigma = np.sqrt(15.0)
x0 = 30.0

#set up initial conditions
y = analytic(x,x0,sigma,v,t)

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x,y,'x-') # numerical data
mpl.plot(x,analytic(x,x0,sigma,v,t),'r-') # analytic data
mpl.show()

yold2 = y
yold  = y
ntmax = 2000
for it in range(ntmax):
    t += dt
    # save previous and previous previous data
    yold2 = yold
    yold  = y

    # get new data; ideally just call a function
    y = yold  + dt*rhs(x,y,v,derivative="central") # FTCS
    #y = yold  + dt*rhs(x,y,v,derivative="left")    # Upwind
    #y = yold  + dt*rhs(x,y,v,derivative="right")   # Downwind
    #y[1:-1] = (yold[:-2]+yold[2:])/2. + dt*rhs(x,y,v,derivative="central")[1:-1] # Lax-Friedrich
    #y = yold2 + 2*dt*rhs(x,y,v,derivative="central") # Leapfrog

    # after update, apply boundary conditions
    y = apply_bcs(x,y) 

    # get analytic result for time t
    yana = analytic(x,x0,sigma,v,t)
    # compute error estimage
    err = np.sum(np.abs(y-yana)/np.abs(yana))
    print "it = ",it,err
    mpl.clf()
    # plot numerical result
    mpl.plot(x,y,'x-')
    # plot analytic results
    mpl.plot(x,yana,'r-')
    mpl.draw()


mpl.show()


