"""
Script to convert the provided spectra from ashely into the same format as the rest of the spectra.
"""

import numpy as np
from astropy.io import fits
# BP Necessary imports.

path = '/d/tel2/brock/Data/Extinction/'
# BP Adding file path to spectra taken by Ashley.

comb_names = ['pk04_a', 'pk04_b', 'pk04_c']
# PB Defining list of files to combine.

comb = 0
# BP Defining empty variable to add combine data to.

# BP Looping over files of the same star to combine into one csv file.
for j in range(len(comb_names)):
	inname = path + comb_names[j] + '_ashley.fits'
	# BP Defining full path to files.
	hdul = fits.open(inname)
	# BP Opening fits file.
	comb += hdul[0].data
	# BP Combining data into one array.

hdul[0].data = comb
# BP Setting data of the last fits file to the newly combined array.
hdul.writeto(path + 'pk04_ashley.fits', overwrite=True)
# BP Write data to new fits file.
# BP Files that were combined have since been deleted for succinctness. 

names = ['pk04', 'pk13', 'pk24', 'pk26', 'pk27']
# BP Defining target names.

for i in range(len(names)):
	inname = path + names[i] + '_ashley_norm.fits'
	# BP Full filepath to read in.
	hdul = fits.open(inname)
	# BP Opening hdu list of defined fits files.
	flux = hdul[0].data
	# BP reading in data.
	hdr = hdul[0].header
	# BP Defining header.
	crval = hdr['CRVAL1']
	crpix = hdr['CRPIX1']
	cdelt = hdr['CD1_1']
	# BP Defining parameters to create wavelength array.
	
	pix = np.linspace(1, len(flux), len(flux))
	# BP Creating array of pixels
	wavelength = crval + cdelt*(pix - crpix)
	# BP Creating array of wavelengths.
	
	flux_err = flux * 0.02
	# BP Using 2 percent error on each data point as an upper limit on Poisson noise.
	
	outdata = np.column_stack((wavelength, flux, flux_err))
	# BP Combining data array as three columns to save to file.
	
	np.savetxt(path + names[i] + '_opt_norm.dat', outdata, header = '     Wavelength                    Flux                  Error')
	# BP Saving file with spectrum information.
