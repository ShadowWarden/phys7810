#!/usr/bin/env python2.7
#
# crank_nicolson.py
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
# <Modified from harmosceuler.py>
#
#
# Reference material at http://docs.scipy.org/doc/numpy/reference/
#

# Import basics to get started
import os, sys, getopt
import math
import string
import numpy
import matplotlib
matplotlib.use("Agg")
matplotlib.rc('axes', hold = True)
from matplotlib import pylab, ticker

flag = 1 # 0 = second derivative in x
		 # 1 = two ODEs for x and v
Npi = 6
Nloops = 45
N = numpy.linspace(2,47,Nloops, True, 'd')
DT = numpy.zeros(Nloops)
Freq = numpy.zeros(Nloops)
n = 0

for numpts in N[0]:
	tspc = numpy.linspace(0., Npi*math.pi, numpts, True, 'd')
	freq = numpy.linspace(0., Npi/2, Npi/2)
	t = tspc[0]
	dt = tspc[1]
	DT[n] = dt

	print "t = ", t
	print "dt = ", dt
	x = numpy.cos(t)
	print x

	fontsz = 20

	om2 = 1.
	pdt2om2 = 1./(1+0.25*dt**2*om2)
	x[0] = 1.
	v = numpy.zeros(int(numpts), 'd')
	print v

	# CN core loop
	for i in range(1,int(numpts)):
		v[i] = (v[i-1]*(1-0.25*om2*dt**2)-om2*dt*x[i-1])*pdt2om2
		x[i] = x[i-1] + 0.5*dt*(v[i] + v[i-1])

	flag = 0
	time = 0
	# Time crossing loop
	for i in range(1,int(numpts)):
		if((x[i] > 0 and x[i-1] < 0) or (x[i] < 0 and x[i-1]>0)):
		# Compute the time delay
			if(((flag%2) == 0) and (flag != 0)):
				timeold = time
				time = tspc[0][i-1] + (tspc[0][i]-tspc[0][i-1])/(numpy.abs(x[i]-x[i-1]))*(x[i-1])
				print(i, time)
				freq[flag/2] = 1./(time-timeold)

			flag += 1
	print(freq)
	freq = freq[1:]
	print(freq.mean())
	Freq[n] = freq.mean()
	n+=1

print(Freq)
Freq_error = abs(Freq - om2/(2*math.pi))
pylab.plot(numpy.log(DT),numpy.log(Freq_error))
# Least squares fit of loglog plot
xdata = numpy.log(DT)
ydata = numpy.log(Freq_error)
A = numpy.vstack([xdata, numpy.ones(len(xdata))]).T
m, c = numpy.linalg.lstsq(A, ydata)[0]
print(m,c)
pylab.plot(xdata,m*xdata+c)
pylab.grid()
pylab.savefig("freqplot.png")

pylab.figure(2)
pylab.plot(t, numpy.cos(t),'b')
pylab.plot(t, x, 'r')
pylab.grid()
pylab.savefig("xho.png")
pylab.close()
