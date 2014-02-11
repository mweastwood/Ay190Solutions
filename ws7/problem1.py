#!/usr/bin/env python

from __future__ import division
from pylab import *

seed(12345) # Seed the random number generator

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
        x = rand(n,2) # Generate the random numbers
        inside_circle = sum((x[:,0]**2 + x[:,1]**2) < 1)
        trial_estimates[j] = 4*inside_circle/n
    estimate[i] = mean(trial_estimates)
    error[i]    = std(trial_estimates)/sqrt(N[i])

figure()
errorbar(log10(N),estimate,yerr=error,fmt='ko')
axhline(y=pi,color='r')
xlabel(r'$\log_{10}N$',fontsize=20)
ylabel(r'$\pi$',fontsize=20)
xlim(0,7)
ylim(pi-2,pi+2)
savefig('problem1.pdf')
close()

