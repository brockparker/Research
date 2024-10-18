import numpy as np
import matplotlib.pyplot as plt

# load a photometric dataset for xxxx
ptimes, mags, merrs = np.loadtxt('/d/users/brock/TESS/TOI_test.csv', skiprows=1,unpack=True,delimiter=',',usecols=[0,3,4])

#ptimes, mags, merrs = np.loadtxt('velcurve.dat', skiprows=1,unpack=True,usecols=[0,1,2])

# look at light curve power
from astropy.stats import LombScargle
frequency = np.linspace(0.035, .05, 10000)
maxp=0
for i in range(1,4):
   power = LombScargle(ptimes,mags,merrs,nterms=i).power(frequency)
   if (np.max(power) > maxp):
      maxp=np.max(power)
      maxind=np.where(power==maxp )
      maxi=int(maxind[0])
      maxfreq=frequency[maxi]
#      print(i,maxi,maxp,maxfreq)
      maxcomp=i
      bestpower=power # power array to plot
print('Max power of {0:5.2f} at freq {1:2.5f} with period {2:2.5f} d with Fourier components= {3:}'.format(maxp,maxfreq,1./maxfreq,maxcomp) )
plt.plot(frequency, bestpower)       
plt.axis([frequency[0],frequency[-1],0,np.max(bestpower)])
plt.show()
