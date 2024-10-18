'''Script to smooth input spectrum using boxcar averaging.'''

import numpy as np
# BP Necessary imports.

path = '/d/tel2/brock/Data/HST/cyc29/data/'
# BP Adding file path to data files.
infiles = ['pk04', 'pk13', 'pk24', 'pk26', 'pk27', 'zeta_oph']
# BP Defining input filenames.

smooth_size = [100, 100, 100, 100, 100, 50]
# BP Defining number of pixels for averaging, approximately 150 AA per block.

for i in range(len(infiles)):
	wavelength, flux, flux_err = np.loadtxt(path + '/' + infiles[i] + '.dat', unpack = True)
	# BP Load in data from each file.
	
	offset = 20
	# BP Defining offset from either end of spectrum to guard against bad data.
	new_size = int(len(wavelength) / smooth_size[i])
	# BP Defining size of new binned data array.
	
	new_wavelength, new_flux, new_flux_err = np.zeros(new_size), np.zeros(new_size), np.zeros(new_size)
	# BP Creating blank arrays to save new data to.
	
	old_flux_err = np.zeros(new_size)
	
	for k in range(new_size):
		new_wavelength[k] = np.mean(wavelength[k*smooth_size[i] + offset:(k+1)*smooth_size[i] + offset])
		# BP Average wavelengths to find center bin wavelengths.
		new_flux[k] = np.mean(flux[k*smooth_size[i] + offset:(k+1)*smooth_size[i] + offset])
		# BP Average fluxes to find center bin fluxes.
		
		### old_flux_err[k] = np.mean(flux_err[k*smooth_size[i] + offset:(k+1)*smooth_size[i] + offset])
		
		new_flux_err[k] = np.sqrt(np.sum(np.square(flux_err[k*smooth_size[i] + offset:(k+1)*smooth_size[i] + offset])) / len(flux_err[k*smooth_size[i] + offset:(k+1)*smooth_size[i] + offset]))
		# BP Combine errors in quadrature and divide error by root N.
		
	outname = path + 'sm' + infiles[i] + '.dat'
	# BP Creating outfile names to save to.

	outdata = np.column_stack((new_wavelength, new_flux, new_flux_err))
	# BP Combining data arrays as three columns to save to file.

	np.savetxt(outname, outdata, header = '      Wavelength                  Flux                  Error')
	# BP Saving file with spectrum information.
