#!/usr/bin/env python

from __future__ import division
from pylab import *
from scipy import integrate
from scipy import interpolate

G = 6.67e-8 # cgs units

# (a) Load and plot the MESA output
r,rho = genfromtxt('presupernova.dat',usecols=(2,4),unpack=True)
rho = rho[r < 1e9]
r   = r  [r < 1e9]

figure()
plot(log10(r),log10(rho),'k-')
xlabel(r'$\log_{10}r/{\rm cm}$',fontsize=20)
ylabel(r'$\log_{10}\rho/{\rm g}\,{\rm cm}^{-3}$',fontsize=20)
savefig('density.png')
close()

# (b) Interpolate onto a grid
f   = interpolate.InterpolatedUnivariateSpline(r,rho)
r   = linspace(r.min(),r.max(),100)
rho = f(r)

# (c) Solve for the potential using an ODE method
def ode_method(r,rho):
    f = interpolate.InterpolatedUnivariateSpline(r,rho)

    def RK4(x,r,dr,ode):
        k1 = dr*ode(r,     x)
        k2 = dr*ode(r+dr/2,x+k1/2)
        k3 = dr*ode(r+dr/2,x+k2/2)
        k4 = dr*ode(r+dr,  x+k3)
        return x + (k1+2*k2+2*k3+k4)/6

    # Define the ODE (dM/dr, d2phi/dr2, dphi/dr)
    dr  = r[2]-r[1]
    ode = lambda r,x: array([4*pi*f(r)*r**2, \
                             4*pi*G*f(r)-2/r*x[1], \
                             x[1]])

    # Integrate!
    x = zeros((len(r),3))
    x[0,:] = (0,0,0) # Initial conditions (M,dphi/dr,phi)
    for i in range(len(r)-1):
        x[i+1,:] = RK4(x[i,:],r[i],dr,ode)
    M   = x[-1,0]
    phi = x[ :,2]

    # Match the boundary conditions
    phi -= G*M/r.max() + phi[-1]
    return phi

# - test case (constant density)
figure()
test_r   = linspace(1e-10,1)
test_rho = ones(test_r.shape)
phi = ode_method(test_r,test_rho)
plot(test_r,2/3*pi*G*(test_r**2-3),'k-',label='analytic solution')
plot(test_r,phi,'ko',label='numerical solution')
xlabel(r'$r$',fontsize=20)
ylabel(r'$\phi$',fontsize=20)
legend(loc='best')
savefig('test_potential.png')
close()

# - presupernova case
figure()
phi = ode_method(r,rho)
plot(r,phi,'ko',label='numerical solution')
xlabel(r'$r / cm$',fontsize=20)
ylabel(r'$\phi / {\rm cgs}\,{\rm units}$',fontsize=20)
savefig('presupernova_potential.png')
close()

