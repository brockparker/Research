'''Script to read in HST STIS spectrum fits files and output wavelength, flux, and flux error in a data file.'''

import numpy as np
from astropy.io import fits 
# BP Necessary imports.

path = '/d/tel2/brock/Data/HST/cyc29/data/'
# BP Adding file path to data fits files.

root = ['oemj01010','oemj02010','oemj03010','oemj04010','oemj05010']
# BP Defining input filenames.
 
outname = ['pk04','pk13','pk24','pk26','pk27'] 
# BP Defining output filenames.

for i in range(len(root)):
	inname = path + '/' + root[i] + '/' + root[i] + '_x1d.fits'
	# BP Full inname for fits files.
	hdul = fits.open(inname)
	# BP Opening hdu list of defined fits file.
	data = hdul[1].data
	# BP Extracting data from proper hdu in list.
	
	wavelengths = np.array(data[0][2])
	# BP Extracting wavelengths.
	flux = np.array(data[0][6])
	# BP Extracting fluxes.
	flux_err = np.array(data[0][7])
	# BP Extracting flux errors.
	
	outdata = np.column_stack((wavelengths, flux, flux_err))
	# BP Combining data arrays as three columns to save to file.

	np.savetxt(path + '/' + outname[i] + '.dat', outdata, header = '      Wavelength                  Flux                  Error')
	# BP Saving file with spectrum information.
  