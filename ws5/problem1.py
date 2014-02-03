#!/usr/bin/env python 
from __future__ import division
from pylab import *
from astropy.io import ascii

def linear_regression(x,y,sigma=None):
    if sigma is None:
        sigma = ones(len(y))

    sumx = sum(x/sigma**2)
    sumy = sum(y/sigma**2)
    sumx2 = sum(x**2/sigma**2)
    sumxy = sum(x*y/sigma**2)
    s = sum(1/sigma**2)

    m = (s*sumxy - sumx*sumy) / \
        (s*sumx2 - sumx**2)
    b = (sumy*sumx2 - sumx*sumxy) / \
        (s*sumx2 - sumx**2)

    return (m,b)

# (a) Read in the data
data = ascii.read('m_sigma_table.dat',readme='m_sigma_ReadMe.dat')  

sigma  = array(data['sigma*'])
logM   = array(data['logM'])
esigma = array(data['e_sigma*'])
elogM  = array(data['e_logM'])

# (b) Compute a linear regression fit (ignoring errors)
x = log10(sigma)
y = logM
m1,b1 = linear_regression(x,y)

# (c) Compute a linear regression fit (including errors)
xerr = esigma/(log(10)*sigma)
yerr = elogM**2
err  = sqrt(yerr**2 + (m1*yerr)**2)
m2,b2 = linear_regression(x,y,err)

# Plot the results
figure(figsize=(6,6))
t = linspace(min(x),max(x))
errorbar(x,y,xerr=xerr,yerr=yerr,fmt='ko',capsize=0)
plot(t,m1*t+b1,'r-',label='without errors')
plot(t,m2*t+b2,'b-',label='including errors')
legend(loc='best')
xlabel(r'$\log_{10}(\sigma_*/{\rm km}\,{\rm s}^{-1})$',fontsize=20)
ylabel(r'$\log_{10}(M_{\rm BH}/M_\odot)$',fontsize=20)
savefig('problem1.pdf')
close()

