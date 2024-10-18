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

name = 'zetaOph'
# BP Defining name for plotting and saving.

wav_snr, flux_snr, err_snr = red_spec.read_spectrum(name, '/d/tel2/brock/Data/Extinction/')
wav_snr, flux_snr = u.Quantity(wav_snr, u.AA), u.Quantity(flux_snr, u.erg / (u.s * u.cm**2 * u.AA))
# BP Reading in full spectrum and converting to proper units.

### sigma = 2
# BP Defining sigma for signal to noise cut-off.

### snr, snr_lims = red_spec.snr(['Zeta Oph'], np.expand_dims(wav_snr, axis = 0), np.expand_dims(flux_snr, axis = 0), sigma)
# BP Calculating lower SNR limits for fitting.

wavelength, flux, flux_err = red_spec.read_spectrum(name.replace('z','sedz'), '/d/tel2/brock/Data/Extinction/')
wavelength, flux, flux_err = u.Quantity(wavelength, u.AA), u.Quantity(flux, u.erg / (u.s * u.cm**2 * u.AA)), u.Quantity(flux_err, u.erg / (u.s * u.cm**2 * u.AA))
# BP Reading in spectrum from specified star.

ind = (wavelength < 30000*u.AA) & (wavelength > 1000*u.AA)
wavelength, flux, flux_err = wavelength[ind], flux[ind], flux_err[ind]
# BP Removing lower and upper wavelengths limits as determined by SNR.

temperature = 35000
gravity = 3.75
metallicity = -0.0
smooth_factor = 100
# BP Defining stellar atmosphere model parameters.

cs_smooth = red_spec.get_tlusty_sed(temperature, gravity, metallicity, smooth_factor)
# BP Reading in interpolated PHOENIX stellar model spectrum.

model_flux = cs_smooth(wavelength)
# BP Evaluating spectrum at observed wavelengths.

model_av, model_rv, model_av_err, model_rv_err, reddened_model_flux, reddened_flux_lower, reddened_flux_upper = red_spec.fit_ext_lsq(wavelength.value, flux.value, flux_err.value, model_flux, red_model)
# BP Fitting reddened stellar spectrum with provided stellar model using least squares method.

chi_squared, red_chi_squared = red_spec.chi_squared(flux.value, reddened_model_flux, flux_err.value)
# BP Retrieving chi squared.


red_spec.plot_best_fit(name.replace('z','Z').replace('aO','a O'), wavelength.value, flux.value, flux_err.value, reddened_model_flux, reddened_flux_lower, reddened_flux_upper, red_chi_squared, temperature, gravity, metallicity, str(red_model.__name__), model_av, model_rv, model_av_err, model_rv_err)
# BP Plotting best fit reddened model plus errors.