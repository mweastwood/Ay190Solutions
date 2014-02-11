#!/usr/bin/env python

from __future__ import division
from pylab import *

seed(12345) # Seed the random number generator

# Define the function and the ranges from which to
# draw random numbers
f = lambda x: x**2 + 1
x_range = (2, 3); dx = x_range[1]-x_range[0]
y_range = (5,10); dy = y_range[1]-y_range[0]

N = logspace(1,6,6) # Total number of random samples
Ntrials = 10        # Number of trials to do with each N

estimate = zeros(len(N))
error    = zeros(len(N))
for i in range(len(N)):
    trial_estimates = zeros(Ntrials)
    for j in range(Ntrials):
        n = round(N[i]/Ntrials) # Note that because we are doing this Ntrials-times,
                                # we need to reduce the number of random numbers
                                # generated (to be honest!)
        x = rand(n)*dx+x_range[0] # Generate the random x-values
        y = rand(n)*dy+y_range[0] # Generate the random y-values
        below_curve = sum(f(x) - y > 0)
        trial_estimates[j] = y_range[0]*dx + below_curve/n*dx*dy
    estimate[i] = mean(trial_estimates)
    error[i]    = std(trial_estimates)/sqrt(N[i])

figure()
errorbar(log10(N),estimate,yerr=error,fmt='ko')
axhline(y=22/3,color='r')
xlabel(r'$\log_{10}N$',fontsize=20)
ylabel(r'$\int_2^3(x^2+1)\,{\rm d}x$',fontsize=20)
xlim(0,7)
ylim(22/3-2,22/3+2)
savefig('problem3.pdf')
close()

