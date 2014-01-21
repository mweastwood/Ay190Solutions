#!/usr/bin/env python

from __future__ import division
from pylab import *

# Do the unstable calculation
four_thirds     = float32(4/3)
thirteen_thirds = float32(13/3)

x = zeros(16,dtype=float32)
x[0] = 1
x[1] = 1/3

for n in arange(1,15):
    x[n+1] = thirteen_thirds*x[n] - four_thirds*x[n-1]

# Calculate what the result should have been
n = arange(16)
y = (1/3)**n

# Plot the absolute error
figure()
plot(n,log10(abs(x-y)),'ko')
xlabel(r'$n$',fontsize=20)
ylabel(r'$\log_{10}|\Delta x_n|$',fontsize=20)
savefig('problem1_absoluteerror.pdf')
close()

# Plot the relative error
figure()
plot(n,log10(abs(x-y)/y),'ko')
xlabel(r'$n$',fontsize=20)
ylabel(r'$\log_{10}|\Delta x_n/x_n|$',fontsize=20)
savefig('problem1_relativeerror.pdf')
close()

