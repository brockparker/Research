'''Plots all input stellar spectra on one stacked plot.'''

import os
os.chdir('/d/users/brock/python/')
import extinction
import numpy as np
import astropy.units as u
from scipy.ndimage import gaussian_filter
# BP Necessary imports.

if __name__ == '__main__':
	stars = ['zetaOph']
	temps = ['34000 K']
	# BP Defining list of star names to import.
	
	red_spec = extinction.reddening()
	# BP Creating fake class to store the spectra.
	
	basedir = '/d/tel2/brock/Data/HST/cyc29/data/'
	# BP Defining base directory where all data is stored.
	
	filt_smooth = 0.5
	# Defining std of Gaussian for smoothing kernel.
	 
	for i in range(len(stars)):
		wavelength, flux, flux_err = red_spec.read_spectrum(stars[i], basedir)
		# BP Reading in fluxes from each data file.

	wavelength, flux, flux_err = np.expand_dims(wavelength, axis=0), np.expand_dims(flux, axis=0), np.expand_dims(flux_err, axis=0)

	wavelengths = u.Quantity(wavelength, u.AA)
	fluxes = u.Quantity(flux, u.erg / (u.s * u.cm**2 * u.AA))
	flux_errs = u.Quantity(flux_err, u.erg / (u.s * u.cm**2 * u.AA))
	# BP Converting quanitities into proper units.	

	fluxes[fluxes < 0] = 0
	#fluxes = gaussian_filter(fluxes, filt_smooth) * u.erg / (u.s * u.cm**2 * u.AA)
	# BP Setting all negative fluxes to 0 and convolving with a Gaussian of sigma 1.

	lower, upper = 1400*u.AA, 3800*u.AA
	# BP Defining upper and lower limits for plotting.
	out_stars = [a.replace('z','Z').replace('aO','a O') for a in stars]
	# BP Removing 2021 from strings for ease of reading.
	
	lines =  {'Mg II' :   	[2799, 0, 0]*u.AA,
			  'Fe III':	[2555, 0, 0]*u.AA,
			  ' Fe III ':	[2905, 100, 0.1]*u.AA,
			  '   Fe III Forest':	[1854, +210, 0.05]*u.AA,
			  '':		[2100, -35, 0.05]*u.AA}
  	# BP Defining line names, locations, x offsets, and y offsets for plotting.
	
	red_spec.plot_raw(wavelengths, fluxes, flux_errs, out_stars, lower, upper, temperatures = temps, lines = lines, save_file = '/d/users/brock/paper_figures/spectrum_smooth_zeta.pdf')
	# BP Plotting raw spectra.
