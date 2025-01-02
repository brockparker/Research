import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from astropy.io import fits
import csv

t = array([0.000]*10000)

for x in range(0, 10000):
	t[x] = 0.001*x
	
A1=2
A2=-2
A3=-2
A4=1
A5=1
A6=-1
P1=0
P2=-1
P3=pi/4
P4=-0.5
P5=1.5
P6=6
	
I1=A1*np.sin(((2*np.pi*(t-15))/10)+P1)
I2=A2*np.sin(((2*2*np.pi*(t-15))/10)+P2)
I3=A3*np.sin(((2*3*np.pi*(t-15))/10)+P3)
I4=A4*np.sin(((2*4*np.pi*(t-15))/10)+P4)
I5=A5*np.sin(((2*5*np.pi*(t-15))/10)+P5)
I6=A6*np.sin(((2*6*np.pi*(t-15))/10)+P6)

I=I1+I2+I3+I4+I5+I6


plt.plot(t, I, 'k.')
plt.show()
