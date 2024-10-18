'''Plots all input stellar spectra on one stacked plot.'''

import os
os.chdir('/d/users/brock/python/')
import extinction
import numpy as np
import astropy.units as u
from scipy.ndimage import gaussian_filter
from astropy.visualization import quantity_support
# BP Necessary imports.

if __name__ == '__main__':
	quantity_support()
	
	stars = ['pk04','pk13','pk24','pk26','pk27']
	temps = ['5750 K','6250 K','5750 K','6250 K','6500 K']
	# BP Defining list of star names to import.
	
	red_spec = extinction.reddening()
	# BP Creating fake class to store the spectra.
	
	basedir = '/d/tel2/brock/Data/HST/cyc29/data/'
	# BP Defining base directory where all data is stored.
	
	filt_smooth = 1
	# Defining std of Gaussian for smoothing kernel.
	 
	for i in range(len(stars)):
		if i == 0:
			wavelength, flux, flux_err = red_spec.read_spectrum(stars[i], basedir)
			# BP Reading in fluxes from each data file.
			
		else:
			wave_foo, flux_foo, err_foo = red_spec.read_spectrum(stars[i], basedir)
			# BP Reading in fluxes from each data file.
			
			wavelength = np.vstack((wavelength, wave_foo))
			flux = np.vstack((flux, flux_foo))
			flux_err = np.vstack((flux_err, err_foo))
			# BP Stacking spectrum from each star into format expected by extinction.plot_raw().
				
	wavelengths = u.Quantity(wavelength, u.AA)
	fluxes = u.Quantity(flux, u.erg / (u.s * u.cm**2 * u.AA))
	flux_errs = u.Quantity(flux_err, u.erg / (u.s * u.cm**2 * u.AA))
	# BP Converting quanitities into proper units.	
			
	snr, snr_lims = red_spec.snr(stars, wavelengths, fluxes, 2)
	# BP Calculating limits where SNR is less than 2.

	fluxes[fluxes < 0] = 0
	fluxes = gaussian_filter(fluxes, filt_smooth) * u.erg / (u.s * u.cm**2 * u.AA)
	# BP Setting all negative fluxes to 0 and convolving with a Gaussian of sigma 1.

	lower, upper = 1500*u.AA, 3155*u.AA
	# BP Defining upper and lower limits for plotting.
	out_stars = [a.replace('pk','PK-') for a in stars]
	# BP Removing 2021 from strings for ease of reading.
	
	lines =  {'C I' :     [1931, 0, 0]*u.AA,
			  'Fe II/I' : [2745, -200, 0.5]*u.AA,
			  'Mg II' :   [2799, -90, 0.25]*u.AA,
			  'Mg I' :    [2852, 0, 0.25]*u.AA,
			  'Si I' :    [2881, 100, 0.25]*u.AA}
	# BP Defining line names, locations, x offsets, and y offsets for plotting.
	
	red_spec.plot_raw(wavelengths, fluxes, flux_errs, out_stars, lower, upper, temperatures = temps, limits = snr_lims, lines = lines, save_file = '/d/users/brock/paper_figures/spectrum_smooth_all.pdf')
	# BP Plotting raw spectra.
