import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from astropy.io import fits
import csv

hdu = fits.open('TOI_270c.fits')
hdr1 = hdu[1].header
hdr2 = hdu[2].header
data1 = hdu[1].data
data2 = hdu[2].data
data3 = hdu[3].data

time = data2['TIME']
light = data2['LC_DETREND']
err = data2['LC_INIT_ERR']

size = np.size(time)

row = [0]*size
row_list = []

for x in range(0, size):
	row[x] = data2[x]
	if str(row[x][4]) != "nan":
		row_list.append(row[x])
		

with open('/d/users/brock/TESS/TOI_270c.csv', 'w', newline='') as file:
     writer = csv.writer(file)
     writer.writerows(row_list)

ptimes, mags, merrs = np.loadtxt('/d/users/brock/TESS/TOI_1581_ASASSN.csv', skiprows=1,unpack=True,delimiter=',',usecols=[0,3,4])
period = 24.44803
phase = (np.mod((ptimes-2459000),period))/period

plt.plot(phase, mags, 'b.', linestyle='')
plt.axis([0,1,10.65,10.3])
plt.show()
