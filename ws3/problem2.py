#!/usr/bin/env python

from __future__ import division
from pylab import *
from scipy import special

# constants in cgs units
c    = 2.998e10     # speed of light
kB   = 1.38065e-16  # Boltzmann constant
hbar = 1.054572e-27 # reduced Planck constant

MeV_to_ergs = 1.602177e-6 # factor to convert units of energy from MeV to ergs
T = 20*MeV_to_ergs/kB     # temperature (kT = 20MeV)
integrand = lambda x: x**2/(exp(x)+1)

# (a) Use Gauss-Laguerre Quadrature to determine the total
#     number density of electrons
[laguerre_roots,laguerre_weights] = special.l_roots(20)
weight   = lambda x: exp(-x)
integral = sum(laguerre_weights*integrand(laguerre_roots)/weight(laguerre_roots))

n_electrons = 8*pi*(kB*T)**3/(2*pi*hbar*c)**3 * integral
print 'The number density of electrons is %.4e/cm^3' % n_electrons

# (b) Use Gauss-Legendre Quadrature to determine the spectral
#     distribution of the electrons
[legendre_roots,legendre_weights] = special.p_roots(20)

dE = 5; Emax = 200
E_bins = linspace(0,Emax,Emax/dE+1) # MeV
x_bins = E_bins*MeV_to_ergs/(kB*T)
n_bins = zeros(len(E_bins)-1) # number density of electrons in each energy bin

for i in range(len(E_bins)-1):
    a = x_bins[i]; b = x_bins[i+1]
    transformed_integrand = lambda x: (b-a)/2*integrand((b-a)/2*x+(a+b)/2)
    integral = sum(legendre_weights*transformed_integrand(legendre_roots))
    n_bins[i] = 8*pi*(kB*T)**3/(2*pi*hbar*c)**3 * integral / dE
print 'The integral of the spectral distribution is %.4e/cm^3' % sum(n_bins*dE)

figure()
plot((E_bins[1:]+E_bins[:-1])/2,n_bins,'k-')
xlabel(r'$E / {\rm MeV}$',fontsize=20)
ylabel(r'${\rm d}n_{e^\pm}/{\rm d}E$',fontsize=20)
savefig('problem2.pdf')
close()

