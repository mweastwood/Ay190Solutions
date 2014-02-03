#!/usr/bin/env python

from __future__ import division
from pylab import *
from timeit import timeit

import student_implementations
students = array(['anthony','cutter','danielk','daniel','david','donal','gina', \
                  'hannah','john','matthias','mee','sarah','scott','xiangcheng'])

# Check DFT correctness
x = randn(50)
y = fft(x)
for student in students:
    module = vars(student_implementations)[student]
    z = module.dft(x)
    if sum(abs(y-z)) > 1e-10:
        print '%s\'s dft is incorrect.' % student

# Time each DFT
N = arange(10,100,10)
times = [[timeit('dft(pylab.randn(%d))' % n, number = 100, \
                 setup='import pylab; from student_implementations.%s import dft' % student) \
          for n in N] for student in students]
times = array(times)

# Find the fastest algorithms
ind = argsort(times[:,-1])
print 'The fastest algorithms, in order are:'
print students[ind]

# Plot the results
figure()
colors = cm.hsv(linspace(0,1,len(students)))
for i,student in enumerate(students):
    plot(log10(N),log10(times[i,:]),
         label=student,color=colors[i,:])
xlabel(r'$\log_{10}N$',fontsize=20)
ylabel(r'$\log_{10}(t/{\rm seconds})$',fontsize=20)
legend(loc='best')
savefig('timing.pdf')
close()

