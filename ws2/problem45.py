#!/usr/bin/env python

from __future__ import division
from pylab import *
from scipy import interpolate

def construct_lagrange_interpolant(xdata,ydata):
    '''
    This function returns another function that calculates
    the Lagrange interpolant of (xdata,ydata)
    '''
    def interpolant(x):
        y = zeros(len(x))
        for xi,yi in zip(xdata,ydata):
            L = ones(len(x))
            for xj, yj in zip(xdata,ydata):
                if xj == xi and yj == yi: continue
                L *= (x-xj)/(xi-xj)
            y += L*yi
        return y
    return interpolant

def construct_linear_interpolant(xdata,ydata):
    '''
    This function returns another function that calculates
    the piecewise linear interpolant of (xdata,ydata)
    '''
    def interpolant(x):
        y = zeros(len(x))
        for i,xi in enumerate(x):
            # Find the nearest points on either side of xi
            xj = xdata[xdata <= xi].max()
            xk = xdata[xdata >= xi].min()
            yj = ydata[xdata == xj]
            yk = ydata[xdata == xk]

            # Linearly interpolate to xi
            if xj == xk: # xi == xj == xk
                y[i] = yj
            else:
                y[i] = yj + (yk-yj)/(xk-xj)*(xi-xj)
        return y
    return interpolant

def construct_quadratic_interpolant(xdata,ydata):
    '''
    This function returns another function that calculates
    the piecewise quadratic interpolant of (xdata,ydata)
    '''
    def interpolant(x):
        y = zeros(len(x))
        for i,xi in enumerate(x):
            # Find the nearest points on either side of xi
            xlo = sort(xdata[xdata <= xi])[::-1]
            xhi = sort(xdata[xdata >= xi])

            if len(xhi) < 2: # Take two points from xlo and one from xhi
                xj,xk = xlo[:2]
                xl    = xhi[0]
            else:            # Take one point from xlo and two from xhi
                xj    = xlo[0]
                xk,xl = xhi[:2]
            yj = ydata[xdata == xj]
            yk = ydata[xdata == xk]
            yl = ydata[xdata == xl]

            # Interpolate to xi
            if   xi == xj:
                y[i] = yj
            elif xi == xk:
                y[i] = yk
            elif xi == xl:
                y[i] = yl
            else:
                y[i] = (xi-xk)*(xi-xl)/(xj-xk)/(xj-xl)*yj \
                     + (xi-xj)*(xi-xl)/(xk-xj)/(xk-xl)*yk \
                     + (xi-xj)*(xi-xk)/(xl-xj)/(xl-xk)*yl
        return y
    return interpolant

def construct_cubic_hermite_interpolant(xdata,ydata):
    '''
    This function returns another function that calculates
    the cubic Hermite interpolant of (xdata,ydata)
    '''
    def interpolant(x):
        y = zeros(len(x))
        phi0 = lambda z: 2*z**3 - 3*z**2 + 1
        phi1 = lambda z:   z**3 - 2*z**2 + z

        # Estimate the derivative at each point with central differencing
        deriv = (ydata[:-2]-ydata[2:])/(xdata[:-2]-xdata[2:])
        xdata_minus_endpoints = xdata[1:-1]
        ydata_minus_endpoints = ydata[1:-1]

        if (x < xdata_minus_endpoints.min()).any() or (x > xdata_minus_endpoints.max()).any():
            raise ValueError('x must be within the interpolation region')

        for i,xi in enumerate(x):
            # Find the nearest points on either side of xi
            xj = xdata_minus_endpoints[xdata_minus_endpoints <= xi].max()
            xk = xdata_minus_endpoints[xdata_minus_endpoints >= xi].min()
            yj = ydata_minus_endpoints[xdata_minus_endpoints == xj]
            yk = ydata_minus_endpoints[xdata_minus_endpoints == xk]
            zj = deriv[xdata_minus_endpoints == xj]
            zk = deriv[xdata_minus_endpoints == xk]

            # Interpolate to xi
            if xj == xk: # xi == xj == xk
                y[i] = yj
            else:
                z = (xi-xj)/(xk-xj)
                y[i] = yj*phi0(z) + yk*phi0(1-z) + zj*(xk-xj)*phi1(z) - zk*(xk-xj)*phi1(1-z)
        return y
    return interpolant

# Load the data
t,mag = genfromtxt('cepheid_lightcurve.txt',unpack=True)

# Construct all the interpolants
lagrange_interpolant  =      construct_lagrange_interpolant(t,mag)
linear_interpolant    =        construct_linear_interpolant(t,mag)
quadratic_interpolant =     construct_quadratic_interpolant(t,mag)
hermite_interpolant   = construct_cubic_hermite_interpolant(t,mag)
spline_interpolant    = interpolate.InterpolatedUnivariateSpline(t,mag,k=3)

# Plot!
figure(figsize=(12,8))
t1 = linspace(0,1,1e3)
t2 = linspace(0.2,0.8,1e3)

colors = cm.jet_r(linspace(0,1,5))
plot(t,mag,'ko')
plot(t1, lagrange_interpolant(t1),color=colors[0],label='Lagrange')
plot(t1,   linear_interpolant(t1),color=colors[1],label='piecewise linear')
plot(t1,quadratic_interpolant(t1),color=colors[2],label='piecewise quadratic')
plot(t2,  hermite_interpolant(t2),color=colors[3],label='piecewise cubic Hermite')
plot(t1,   spline_interpolant(t1),color=colors[4],label='cubic spline')
xlabel(r'$t / {\rm days}$',fontsize=20)
ylabel(r'$m$',fontsize=20)
leg = legend(loc='best')
leg.get_frame().set_alpha(0.5)
xlim(0,1); ylim(0,0.8)
gca().invert_yaxis()
savefig('problem45.pdf')
close()

