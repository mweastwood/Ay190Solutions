#!/usr/bin/env python

from __future__ import division
import numpy as np
import scipy as sp


# global constants
ggrav = 6.67e-8
msun  = 1.99e33

# EOS parameters
# for white dwarfs:
polyG = 4.0/3.0
polyK = 1.244e15*0.5**polyG


#######################################
# function definitions
def tov_RHS(rad,p,m):
    
    # RHS function

    # invert EOS to get density (needed in rhs)
    rho = (p/polyK)**(1.0/polyG)

    rhs = np.zeros(2)
    if(rad > 1.0e-10):
        rhs[0] = - ggrav * m / rad**2 * rho
        rhs[1] = 4.0 * np.pi * rad**2 * rho
    else:
        rhs[0] = 0.0
        rhs[1] = 0.0

    return rhs

def tov_integrate_FE(rad,dr,p,m):

    # Forward-Euler Integrator

    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = p
    old[1] = m

    # forward Euler integrator
    # for this, must call RHS routine
    new = old + dr * tov_RHS(rad,p,m)
    
    # assign outputs
    pnew = new[0]
    mnew = new[1]
    
    return (pnew,mnew)

def tov_integrate_RK2(r,dr,p,m):

    dp1,dm1 = dr*tov_RHS(r,     p,      m)
    dp2,dm2 = dr*tov_RHS(r+dr/2,p+dp1/2,m+dm1/2)

    pnew = p+dp2
    mnew = m+dm2

    return (pnew,mnew)

def tov_integrate_RK3(r,dr,p,m):

    dp1,dm1 = dr*tov_RHS(r,     p,          m)
    dp2,dm2 = dr*tov_RHS(r+dr/2,p+dp1/2,    m+dm1/2)
    dp3,dm3 = dr*tov_RHS(r+dr,  p-dp1+2*dp2,m-dm1+2*dm2)

    pnew = p+(dp1+4*dp2+dp3)/6
    mnew = m+(dm1+4*dm2+dm3)/6

    return (pnew,mnew)

def tov_integrate_RK4(r,dr,p,m):

    dp1,dm1 = dr*tov_RHS(r,     p,      m)
    dp2,dm2 = dr*tov_RHS(r+dr/2,p+dp1/2,m+dm1/2)
    dp3,dm3 = dr*tov_RHS(r+dr/2,p+dp2/2,m+dm2/2)
    dp4,dm4 = dr*tov_RHS(r+dr,  p+dp3,  m+dm3)

    pnew = p+(dp1+2*dp2+2*dp3+dp4)/6
    mnew = m+(dm1+2*dm2+2*dm3+dm4)/6

    return (pnew,mnew)

#######################################

# set up grid
npoints = 10000
radmax = 2.0e8 # 2000 km
radius = np.linspace(0.0,radmax,npoints)
dr = radius[1]-radius[0]

# set up variables
press = np.zeros(npoints)
rho   = np.zeros(npoints)
mass  = np.zeros(npoints)

# set up central values
rho[0]   = 1.0e10
press[0] = polyK * rho[0]**polyG
mass[0]  = 0.0

# set up termination criterion
press_min = 1.0e-10 * press[0]

nsurf = 0
for n in range(npoints-1):
    
    (press[n+1],mass[n+1]) = tov_integrate_FE (radius[n],dr,
                                               press[n],mass[n])
    # check for termination criterion
    if(press[n+1] < press_min and nsurf==0):
        nsurf = n

    if(n+1 > nsurf and nsurf > 0):
        press[n+1] = press[nsurf]
        rho[n+1]   = rho[nsurf]
        mass[n+1]  = mass[nsurf]

    # invert the EOS to get density
    rho[n+1] = (press[n+1]/polyK)**(1.0/polyG)


print radius[nsurf]/1.0e5
print mass[nsurf]/msun



