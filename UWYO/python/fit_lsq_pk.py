import os
os.chdir('/d/users/brock/python/')
import extinction
import numpy as np
import astropy.units as u

from dust_extinction.parameter_averages import F19, CCM89, VCG04, G23, O94, F99, D22, F04
# BP Necessary imports.

red_model = CCM89
# BP Defining reddening model to use.

red_spec = extinction.reddening()
# BP Initializing reddening class.

name = 'pk04'
# BP Defining name for plotting and saving.

wav_snr, flux_snr, err_snr = red_spec.read_spectrum(name, '/d/tel2/brock/Data/Extinction/')
wav_snr, flux_snr = u.Quantity(wav_snr, u.AA), u.Quantity(flux_snr, u.erg / (u.s * u.cm**2 * u.AA))
# BP Reading in full spectrum and converting to proper units.

sigma = 2.00
# BP Defining sigma for signal to noise cut-off.

snr, snr_lims = red_spec.snr([name], np.expand_dims(wav_snr, axis = 0), np.expand_dims(flux_snr, axis = 0), sigma)
# BP Calculating lower SNR limits for fitting.

wavelength, flux, flux_err = red_spec.read_spectrum(name.replace('pk','sedPK'), '/d/tel2/brock/Data/Extinction/')
wavelength, flux, flux_err = u.Quantity(wavelength, u.AA), u.Quantity(flux, u.erg / (u.s * u.cm**2 * u.AA)), u.Quantity(flux_err, u.erg / (u.s * u.cm**2 * u.AA))
# BP Reading in spectrum from specified star.

ind = (wavelength > snr_lims[0]) & (wavelength < 30000*u.AA)
wavelength, flux, flux_err = wavelength[ind], flux[ind], flux_err[ind]
# BP Removing lower and upper wavelengths limits as determined by SNR.

temperature_mean = 5750
temperature_sigma = 100
gravity = 4.127
metallicity_mean = -0.218
metallicity_sigma = 0.1
smooth_factor = 400
# BP Defining stellar atmosphere model parameters.

for i in range(1):
	metallicity = metallicity_mean #metallicity_sigma * np.random.randn() + metallicity_mean
	temperature = temperature_mean #temperature_sigma * np.random.randn() + temperature_mean
	
	print(metallicity, temperature)
		
	cs_smooth = red_spec.get_phoenix_sed(temperature, gravity, metallicity, smooth_factor)
	# BP Reading in interpolated PHOENIX stellar model spectrum.
	
	model_flux = cs_smooth(wavelength)
	# BP Evaluating spectrum at observed wavelengths.
	
	model_av, model_rv, model_av_err, model_rv_err, reddened_model_flux, reddened_flux_lower, reddened_flux_upper = red_spec.fit_ext_lsq(wavelength.value, flux.value, flux_err.value, model_flux, red_model)
	# BP Fitting reddened stellar spectrum with provided stellar model using least squares method.
	
	chi_squared, red_chi_squared = red_spec.chi_squared(flux.value, reddened_model_flux, flux_err.value)
	# BP Retrieving chi squared.
	
	outdata = np.column_stack((temperature, gravity, metallicity, chi_squared, model_av, model_rv))
	
	#with open('/d/tel2/brock/Data/Extinction/' + name + '_posterior_' + str(red_model.__name__) + '.txt', 'ab') as f:
	#	np.savetxt(f, outdata)
	
	red_spec.plot_best_fit(name, wavelength.value, flux.value, flux_err.value, reddened_model_flux, reddened_flux_lower, reddened_flux_upper, red_chi_squared, temperature, gravity, metallicity, str(red_model.__name__), model_av, model_rv, model_av_err, model_rv_err)
	# BP Plotting best fit reddened model plus errors.
	
# =============================================================================
# 	best_av = 0.86
# 	best_rv = 1.99
# 	best_av_err = 0.07
# 	best_rv_err = 0.19
# 	
# 	red_spec.plot_best_fit(name, wavelength.value, flux.value, flux_err.value, reddened_model_flux, reddened_flux_lower, reddened_flux_upper, red_chi_squared, temperature, gravity, metallicity, str(red_model.__name__), best_av, best_rv, best_av_err, best_rv_err)
# 	# BP Plotting best fit reddened model plus errors.
# =============================================================================
	