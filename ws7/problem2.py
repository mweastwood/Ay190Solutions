#!/usr/bin/env python

from __future__ import division
from pylab import *

def is_there_a_duplicate(l):
    if len(unique(l)) < len(l):
        return True
    else:
        return False

def probability_of_same_birthday(N):
    if N >= 365: return 1
    not_probability = 1
    for i in range(1,N):
        not_probability *= (365-i)/365
    return 1-not_probability

seed(12345)      # Seed the random number generator
N = arange(2,30) # Total number of people

# Calculate the probability 2 people share a birthday
# using Monte Carlo methods
trials = 10000 # Number of trials to run for each number of people
fraction = zeros(len(N))

for i,n in enumerate(N):
    count = 0
    for j in range(trials):
        birthdays = randint(0,365,n)
        if is_there_a_duplicate(birthdays):
            count += 1
    fraction[i] = count/trials

# Calculate the probability 2 people share a birthday
# using analytic methods
analytic_fraction = array([probability_of_same_birthday(n) for n in N])

figure()
plot(N,fraction,'k-',label='Monte Carlo')
plot(N,analytic_fraction,'r--',label='analytic')
legend(loc='best')
xlabel('number of people')
ylabel('probability of at least 2 sharing a birthday')
savefig('problem2.pdf')
close()

