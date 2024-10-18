import os
import sys
import fnmatch
import itertools
import numpy as np
from astropy import stats
from astropy.io import fits
from tkinter import Tk
from tkinter.filedialog import askdirectory
from collections import namedtuple

###########################################################################
# Input Target Data Directory When Prompted

# Indicate Naming Convention for Darks, Biases, Flats, and Science Images When Prompted

# This script takes input data, including calibration images, and automatically reduces data to a point acceptable for light curve fitting
###########################################################################

###########################################################################
# Defined Functions Used in script


def find(pattern, loc):
	'''Locates files with a given pattern and outputs them as a list of strings'''
	result = []
	for root, dirs, files in os.walk(loc):
		for name in files:
			if fnmatch.fnmatch(name, pattern):
				result.append(name)
	return result


def sig_clip(array, sigma):
	'''Sigma clips images pixel by pixel and outputs 3D array with nans in place of clipped pixels'''
	iterations = 0
	clipped = 1
	pixels = 11
	array = np.array(array, dtype=float)
	while(pixels > 10):
		rms = np.nanstd(array, axis = 2)
		avg = np.nanmean(array, axis = 2)
		med = np.nanmedian(array, axis = 2)
		clipped = 0
		iterations += 1
		for i, j, k in itertools.product(range(array.shape[0]), range(array.shape[1]), range(array.shape[2])):
			if array[i,j,k] > med[i,j]+sigma*rms[i,j]:
				array[i,j,k] = np.nan
				clipped += 1
			elif array[i,j,k] < med[i,j]-sigma*rms[i,j]:
				array[i,j,k] = np.nan
				clipped += 1
		print("Clipped Pixels:" + str(clipped))
		pixels = clipped
	print("Number of Iterations:" + str(iterations))
	return array


###########################################################################

''''''''''''''''''''''''''''''
'''Bias Processing'''
''''''''''''''''''''''''''''''

path = input('Input Data Directory:')
# prompts the user to input the directory where the data is stored

print('Processing Biases...')

bias_in = find('b-*.fit', path)
# finds all biases, assumming common naming convention b-0???.fit

if(len(bias_in) == 0):
	print('No biases found. Recheck directory and try again. Exiting program.')
	sys.exit()
# exits program if there are no biases in the list

for i in range(len(bias_in)):
	foo1 = fits.open(path + '/' + bias_in[i])
	if (i == 0):
		foo2 = foo1[0].data.reshape((2048,2048,-1))
		biashdr = foo1[0].header
		continue
	foo2 = np.concatenate((foo2, foo1[0].data.reshape((2048,2048,-1))), axis=2)
	foo1.close()
# combine all biases into 3d array for averaging

foo3 = sig_clip(foo2, sigma = 3)
# sigma clip 3d bias array
#foo3 = foo2

'''''''''
Median with rejection

My sheet says average with rejection
'''''''''
masterbias = np.nanmedian(foo3, axis = 2)
# average all of biases using median into masterbias

fits.writeto(path + '/' + 'masterbias.fits', masterbias, biashdr, overwrite=True)
# writes array of averaged biases to fits file

bias_tosubtract = find('*.fit', path)
# make a list of all files to subtract biases from

for n in bias_tosubtract:
	foo4 = fits.open(path + '/' + n)
	foo5 = np.subtract(foo4[0].data, masterbias)
	fname = n.replace('.fit', '_b.fits')
	fits.writeto(path + '/' + fname, foo5, foo4[0].header, overwrite=True)
	foo4.close()
# subtract biases from each image in input list and rewrite as newimage_b.fits

''''''''''''''''''''''''''''''
'''Dark Processing'''    
''''''''''''''''''''''''''''''

print('Processing Darks...')

dark_in = find('d*_b.fits', path)
# finds all darks, assuming naming convention d-0???.fit

if(len(dark_in) == 0):
	print('No darks found. Recheck directory and try again. Exiting program.')
	sys.exit()
# exits program if there are no darks in the list

for i in range(len(dark_in)):
	foo6 = fits.open(path + '/' + dark_in[i])
	if (i == 0):
		foo7 = foo6[0].data.reshape((2048,2048,-1))
		darkhdr = foo6[0].header
		continue
	foo7 = np.concatenate((foo7, foo6[0].data.reshape((2048,2048,-1))), axis=2)
	foo6.close()
# combine all darks into 3d array for averaging

foo8 = sig_clip(foo7, sigma = 3)
# sigma clip 3d dark array
#foo8 = foo7

'''''''''
Median with rejection

My sheet says median with rejection
'''''''''
masterdark = np.nanmedian(foo7, axis = 2)
# average all of darks using median into masterdark

fits.writeto(path + '/' + 'masterdark.fits', masterdark, darkhdr, overwrite=True)
# writes array of averaged darks to fits file

masterdark_norm = masterdark/darkhdr['EXPTIME']
# scale master dark to 1 second exposure time

dark_tosubtract = find('*_b.fits', path)
# make a list of all files to subtract darks from

for n in dark_tosubtract:
	foo8 = fits.open(path + '/' + n)
	foo8hdr = foo8[0].header
	if (float(foo8hdr['EXPTIME']) == 0.0):
		pass
	else:
		foo9 = np.subtract(foo8[0].data, masterdark_norm*foo8hdr['EXPTIME'])
		fname = n.replace('_b.fits', '_bd.fits')
		fits.writeto(path + '/' + fname, foo9, foo8hdr, overwrite=True)
		foo8.close()
# sacle masterdark to each image expsure time, divide masterdark from each image and rewrite as newimage_bd.fits

''''''''''''''''''''''''''''''
'''Flats Processing'''    
''''''''''''''''''''''''''''''

print('Processing Flats...')

flat_in = find('f*_bd.fits', path)
# finds all flats, assuming naming convention f-0???.fit

if(len(flat_in) == 0):
	print('No flats found. Recheck directory and try again. Exiting program.')
	sys.exit()
# exits program if there are no flats in the list

flat_filt_list = []

class image:
	blank = []

foo11 = image()
foo13 = image()
foo14 = image()
flathdr = image()
masterflat = image()
masterflat_norm = image()
immedian = image()

for i in range(len(flat_in)):							# combine all flats into 3d array for averaging
	foo10 = fits.open(path + '/' + flat_in[i])				# open file
	flat_filt = foo10[0].header['FILTER']					# read in filter of flat file
	
	if (np.average(foo10[0].data) > 50000):					# skip image if its over 50k counts
		continue
	elif (np.average(foo10[0].data) < 4000):
		continue
	#also add if too low
	
	if (flat_filt not in flat_filt_list):					# if there is a new filter		
		flat_filt_list.append(flat_filt)				# add new filter to list of filters
		setattr(foo11, flat_filt, foo10[0].data.reshape((2048,2048,-1)))# set initial array for first image to stack onto
		setattr(flathdr, flat_filt, foo10[0].header)			# set header for master flat to use
		setattr(immedian, flat_filt, np.median(foo10[0].data))		# calculate median for first image	
		# get rid of ^
		continue
	
	setattr(foo11, flat_filt, np.concatenate((getattr(foo11, flat_filt), foo10[0].data.reshape((2048,2048,-1))), axis=2))
	# for subsequent images, concatenate them onto the already existing data to create 3d image array
	if(np.median(foo10[0].data)>getattr(immedian, flat_filt)):		# if medain of current image is larger than median of first image
		setattr(immedian, flat_filt, np.median(foo10[0].data))		# replace filter median for scaling

	foo10.close()

for f in flat_filt_list:
	for i in range(np.size(getattr(foo11, f), axis=2)):
		curmed = np.median(getattr(foo11,f)[:,:,i])
		if (i == 0):
			foo12 = getattr(foo11,f)[:,:,i].reshape((2048,2048,-1))*(getattr(immedian,f)/curmed)
			continue
		foo12 = np.concatenate((foo12, getattr(foo11,f)[:,:,i].reshape((2048,2048,-1))*(getattr(immedian,f)/curmed)), axis = 2)

	setattr(foo13, f, foo12)

	setattr(foo14, f, sig_clip(getattr(foo13, f), sigma = 3))
	#setattr(foo14, f, getattr(foo13, f))	

	setattr(masterflat, f, np.nanmedian(getattr(foo14, f), axis = 2))
	
	fits.writeto(path + '/' + 'masterflat_' + f + '.fits', getattr(masterflat, f), getattr(flathdr, f), overwrite=True) # writes array of averaged darks to fits file

	setattr(masterflat_norm, f, getattr(masterflat, f)/(np.median(getattr(masterflat, f))))
	# scale master flats to unity

flat_tosubtract = find('*_bd.fits', path)
# make a list of all files to divide flats from

for n in flat_tosubtract:
	foo15 = fits.open(path + '/' + n)
	foo15hdr = foo15[0].header
	try:
		f = foo15hdr['FILTER']
	except KeyError:
		continue
	foo16 = np.divide(foo15[0].data, getattr(masterflat_norm, f), out=np.zeros_like(foo15[0].data), where=getattr(masterflat_norm,f)!=0)
	# only divides for pixels where the masterdark is not 0, where it is 0 the output is also 0
	fname = n.replace('_bd.fits', '_bdf.fits')
	fits.writeto(path + '/' + fname, foo16, foo15hdr, overwrite=True)
	foo15.close()

########################
'''Yee Haw everything works lets gooooooooooooo'''
'''just ignore sigma clipping okay its fine'''
########################
