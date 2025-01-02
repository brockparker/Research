import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy.optimize import curve_fit
from astropy.stats import LombScargle
from astropy.io import fits
import csv

#open and define fits file and sections
hdu = fits.open('/d/users/brock/TESS/TESS Fits/1478.fits')
hdr1 = hdu[1].header
hdr2 = hdu[2].header
data1 = hdu[1].data
data2 = hdu[2].data
data3 = hdu[3].data

#name variables for easy access
time = data1['TIME']
light = data1['LC_DETREND']
err = data1['LC_INIT_ERR']

#call size of the data set to append to
size = np.size(light)

#create blank array to append data to
row = [0]*size
row_list = []

#append nonnull values to new array for processing
for x in range(0, size):
	row[x] = data2[x]
	if str(row[x][4]) != "nan":
		row_list.append(row[x])

#write data to new csv for future analysis
with open('/d/users/brock/TESS/TESS Fits/1487.csv', 'w', newline='') as file:
     writer = csv.writer(file)
     writer.writerows(row_list)

# load a photometric dataset from created csv file
ptimes, mags, merrs = np.loadtxt('/d/users/brock/TESS/TESS Fits/1478.csv', skiprows=1,unpack=True,delimiter=',',usecols=[0,3,4])

# look at light curve power
frequency = np.linspace(0.1, 1, 10000)
maxp=0
for i in range(1,8):
   power = LombScargle(ptimes,mags,merrs,nterms=i).power(frequency)
   if (np.max(power) > maxp):
      maxp=np.max(power)
      maxind=np.where(power==maxp )
      maxi=int(maxind[0])
      maxfreq=frequency[maxi]
      maxcomp=i
      bestpower=power # power array to plot
print('Max power of {0:5.2f} at freq {1:2.5f} with period {2:2.5f} d with Fourier components= {3:}'.format(maxp,maxfreq,1./maxfreq,maxcomp) )
plt.plot(frequency, bestpower)       
plt.axis([frequency[0],frequency[-1],0,np.max(bestpower)])
plt.show()

#define period as found by lombscargle
period = 1./maxfreq

#fold the data around the phase
phase = (np.mod((ptimes-2459000),period))/period
phi = 0

#define sine functions to fit data to
def I(phase, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10):
	return A1*np.sin(1*np.pi*(phase-phi))+A2*np.sin(2*np.pi*(phase-phi))+A3*np.sin(3*np.pi*(phase-phi))+A4*np.sin(4*np.pi*(phase-phi))+A5*np.sin(5*np.pi*(phase-phi))+A6*np.sin(6*np.pi*(phase-phi))+A6*np.sin(7*np.pi*(phase-phi))+A8*np.sin(8*np.pi*(phase-phi))+A9*np.sin(9*np.pi*(phase-phi))+A10*np.sin(10*np.pi*(phase-phi))


#fit data to function
popt, pcov = curve_fit(I, phase, mags)


#plot original data and best fit curve
plt.plot(phase, mags, 'b.', markersize = '0.5')
plt.plot(phase, I(phase, *popt), 'r.', label = 'Fit', markersize = '0.5')
plt.tick_params(direction = 'in')
plt.xlabel('Phase')
plt.ylabel('Flux')
plt.title('658 Folded Light Curve ASASSN')
plt.show()



