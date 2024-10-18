"""
Importable class to deredden spectra using different model reddening laws.

Created on Fri Jun 30 22:32:31 2023

@author: brock
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import astropy.units as u
import os
import pickle
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
#from specutils import Spectrum1D
from specutils.fitting import fit_generic_continuum
from specutils.spectra import SpectralRegion
from astropy.visualization import quantity_support
from scipy.interpolate import RegularGridInterpolator
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip
from scipy.signal import find_peaks
# BP Necessary imports.

class reddening(object):
	def __init__(self):
		"""
		 BP Initializes reddening class with proper global variables.
		
		Parameters
		----------
		None.

		Returns
		-------
		None.

		"""
		plt.rc('axes', labelsize=16)
		plt.rc('figure', titlesize=16)
		plt.rc('xtick', labelsize=16)
		plt.rc('ytick', labelsize=16)
		plt.rc('legend', fontsize=16)
		# BP Matploblib stylization parameters.
		
	def read_spectrum(self, name, basedir):
		"""
		Reads in data generated from read_hst_spec in .dat format with three columns of wavelength, flux, and flux error.

		Parameters
		----------
		name : string
			String filename of the star data to read in. Does not include .dat.
		basedir : string
			Directory datafile is located in.

		Returns
		-------
		wavelength : array_like
			Array of wavelengths in AA.
		flux : array_like
			Array of fluxes.
		flux_err : array_like
			Array of flux standard errors.

		"""
		wavelength, flux, flux_err = np.loadtxt(basedir + name + '.dat', unpack=True, skiprows=1)
		# BP Reading in wavelengths, fluxes, and errors from data file.
	
		return wavelength, flux, flux_err
	
	def plot_raw(self, wavelengths, fluxes, flux_errs, stars, lower, upper, temperatures = None, limits = None, lines = None, save_file = None):
		"""
		BP Plots all input raw SEDs with a few defined UV absorption lines.		

		Parameters
		----------
		wavelengths : astropy.units.quantity.Quantity. 
			2D array of input vacuum wavelengths.
		fluxes : astropy.units.quantity.Quantity. 
			2D array of input fluxes.
		flux_errs : astropy.units.quantity.Quantity. 
			2D array of input flux errors.
		stars : list
			List of input stellar names used in plotting.
		lower : int
			Lower wavelength limit to plot.
		upper : int
			Upper wavelength limit to plot.
		temperatures : list, optional
			List of input stellar temperature strings used in plotting text. The default is None.
		limits : boolean, optional
			The SNR limit below which to remove points. The default is False.
		lines : dictionary, optional
			Dictionary of absorption lines, wavelengths, and offsets to plot. The default is None.
		save_file: string, optional
			String containting the location and filename to save file to. The default is None.

		Returns
		-------
		None.

		"""
		fig, ax = plt.subplots(figsize=(6,8), layout='tight')
		# BP Initializing subplot with correct figure size for PK targets.
		# fig, ax = plt.subplots(figsize=(6,4), layout='tight')
		# BP Initializing subplot with correct figure size for Zeta Oph.		
		# BP Note that the given text sizes work best with a horizontal size of 6.
		
		stars = stars[::-1]
		wavelengths = wavelengths[::-1]
		fluxes = fluxes[::-1]
		flux_errs = flux_errs[::-1]
		# BP Inverting input lists to plot higher numbers at the bottom.

		# BP Looping over all stars in inupt list.
		for i in range(len(stars)):
			star = stars[i]
			wavelength = wavelengths[i]
			flux = fluxes[i]
			flux_err = flux_errs[i]
			# BP Extracting stellar parameters.
			
			ind = np.where((wavelength > lower) & (wavelength < upper))
			# BP Creating index to chop off ends of bad data
			
			wavelength, flux, flux_err = wavelength[ind], flux[ind], flux_err[ind]
			# BP Chopping data.

			if limits is not None:
				# BP If limits are set, cut off all points below.
				snr_limit = limits[i]
				# BP Extracting limit from list.
				limit_index = wavelength > snr_limit
				# BP Creating boolean array to index with.
				wavelength, flux, flux_err = wavelength[limit_index], flux[limit_index], flux_err[limit_index]
				# BP Extracting values beyond the wavelength limit.
			
			shift = 1.02
			# BP Defining vertical shift between each plotted spectrum.

			ax.plot(wavelength, np.log10(flux / (u.erg / (u.s * u.cm**2 * u.AA))) + shift*i, linestyle='-', color='b', marker = '')
			# BP Plotting fluxes vs wavelength in log space as blue lines.
			
			if temperatures is not None:
				# BP If temperatures are given, plot them as text.
				temperature = temperatures[::-1][i]
				# BP Extracting proper temperature from list in reverse order.

				ax.text(2400, -15.25+shift*i, star + ',\n$T_{eff}$=' + temperature, size = 18, horizontalalignment='right')
				# BP Writing star name and temperature above each spectrum for PK stars.
				
				# ax.text(2140, -8.825+shift*i, star + ',\n$T_{eff}$=' + temperature, size = 18, horizontalalignment='left')
				# BP Writing star name and temperature above each spectrum for Zeta Oph.
				
			else:
				ax.text(1900, -15+shift*i, star, size = 18, horizontalalignment='right')
				# BP Writing star name above each spectrum.

		if lines is not None:
			# BP Write out each provided line species as dotted lines.
			spectrum = CubicSpline(wavelengths[0], fluxes[0])
			
			min_flux = np.log10(np.min(fluxes[0]).value)
			# BP Defining minimum flux for plotting
						
			# BP Creating cubic spline of lowest spectrum to find y limits for absorption lines.
			for key, item in lines.items():
				value = item[0].value
				offsetx = item[1].value
				offsety = item[2].value
				# BP Extracting dimensionless quantities in terms of units.

				shift_min = 0.55
				# BP Vertical text shift factor for PK stars.
				
				# shift_min = 0.065
				# BP Vertical text shift factor for Zeta Oph.
				
				shift_spline = 0.1
				# BP Line separation factor (from spectrum and text) for PK Stars.
				
				# shift_spline = 0.025
				# BP Line separation factor (from spectrum and text) for Zeta Oph.
				
				ax.plot([value, value], [min_flux + shift_min + offsety, np.log10(spectrum(value)) - shift_spline], color='k', linestyle='-', linewidth = 1)
				ax.plot([value, value + offsetx], [min_flux + shift_min + offsety, min_flux + shift_min + offsety], color = 'k', linewidth = 1, linestyle='-', marker='')
				ax.plot([value + offsetx, value + offsetx], [min_flux + shift_min + offsety, min_flux + shift_min], color='k', linestyle='-', linewidth = 1)
				# BP Plotting small vertical and horizontal lines to indicate absorption lines.
				ax.text(value + offsetx, min_flux + shift_min - 2 * shift_spline, key, fontsize = 14, horizontalalignment = 'center')
				# BP Labelling absorption lines.

		ax.set_xlabel('Wavelength (Ang.)')
		ax.tick_params(axis='both', direction='in', which='both')
		if len(stars) > 1:
			ax.set_ylabel('log F$_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ Ang$^{-1}$) + offset')
		else:
			ax.set_ylabel('log F$_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ Ang$^{-1}$)')
		ax.margins(x = 0, y = 0.015)
		# BP Setting lables and stylization parameters.
		
		if save_file is not None:
			plt.savefig(save_file)
		plt.show()
		# BP Saving and plotting file

	def snr(self, stars, wavelengths, fluxes, value):
		"""
		Quick function to calculate the SNR of all input spectra and return lower limits where the SNR is always greater than value.

		Parameters
		----------
		wavelength : array_like
			Input array of wavelengths in AA.
		flux : array_like
			Input array of fluxes.
		value : array_like
			Input array of flux errors..

		Returns
		-------
		snr : array_like
			2D array of output continuum SNRs.
		snr_lim : array_like
			Array of lower SNR limits.

		"""		
		snr_lims = []
		# BP Creating empty list to store limits in.
		
		for i in range(len(stars)):
			wavelength = wavelengths[i]
			flux = fluxes[i]
			# BP Extracting each star.
			
			run_avg, run_sig = np.array([]), np.array([])
			# BP Defining empty lists to store running averages in.
			
			width = 25
			# BP Defining width in pixels to calculate SNR over.
			
			for k in range(width,len(wavelength)-width):
				run_avg = np.append(run_avg, np.mean(flux[k:width*2+k].value))
				run_sig = np.append(run_sig, np.std(flux[k:width*2+k].value))
				# BP Calculating running average and standard deviation.
			
			if i == 0:
				snr = np.array([run_avg/run_sig])
			else:
				snr = np.vstack((snr, run_avg/run_sig))
			# BP Calculating continuum signal to noise ratio.
			# BP Concatenating all SNRs into a single list.
												
			snr_lims.append(wavelength[width:-width][np.where(snr[i] < value)[0][-1]])
			# Creating list of lower limits where the SNR is greater than 2.
											
		return snr, snr_lims[::-1]
		# BP Reversing limits to match 
	
	def get_phoenix_sed(self, temperature, gravity, metallicity, smooth_factor):
		"""
		Function to return smooted cubic spline PHOENIX model spectrum to evaluate at given wavelengths.

		Parameters
		----------
		temperature : float
			Stellar effective temperature, from 2300K to 12000K.
		gravity : float
			Stellar surface gravity, from 1.0-6.0.
		metallicity : float
			Stellar metallicity, from +0.5 - -2.0.
		smooth_factor : integer
			Number of pixels of gaussian smoothing. Typically 400 for PHOENIX models.

		Returns
		-------
		cs_smooth : scipy.interpolate.CubicSpline
			Smoothed PHOENIX stellar spectrum evaluatable as a cubic spline at any wavelength.

		"""
		pkl_file = '/d/tel2/brock/Catalogs/PHOENIX/lte{:07.15f}-{:.15f}{:+.15f}-{:.0f}.pkl'.format(temperature, gravity, metallicity, smooth_factor)
		# BP Defining name of generated pickle file containing interpolated PHOENIX spectrum cubic spline.
		
		if not os.path.isfile(pkl_file):
			# BP Checking if target parameters have already been interpolated and splined and exist as a pickle file.
			wav_start, wav_end, wav_sep = 1400, 49000, 0.1
			wav_interpolate = np.arange(wav_start, wav_end, wav_sep)	
			# BP Wavelength range to create spectrum spline over with separation.
			
			temps = np.append(np.arange(2300, 7000, 100), np.arange(7000, 12200, 200))
			# BP Creating temperature grid from 2300K to 7000K every 100K and from 7000K to 12000K every 200K.
			metals = np.array([-2.0, -1.0, -0.5, -0.0, +0.5])
			# BP Creating metallicity grid at -2.0, -1.0, -0.5, -0.0, and +0.5 dex.
			gravs = np.arange(1.0, 6.5, 0.5)
			# BP Creating gravity grid from 1.0 to 6.5 every 0.5.
						
			if (temperature not in temps) or (gravity not in gravs) or (metallicity not in metals):
				# BP Checking if the supplied parameters are one of the pregenerated PHOENIX models.
				interp = self.interpolate_phoenix()
				model_flux = interp([temperature, gravity, metallicity]).ravel()
				model_flux = u.Quantity(model_flux, u.erg / (u.s * u.cm**2 * u.cm))
				# BP Interpolating existing PHOENIX models if the desired target parameters are not an existing model.
				
			else:
				model_flux_file = '/d/tel2/brock/Catalogs/PHOENIX/lte{:05.0f}-{:.2f}{:+.1f}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(temperature, gravity, metallicity)
				# BP Specifying file containing spectrum for given stellar parameters.
				
				model_flux, flux_header = fits.getdata(model_flux_file, header=True)
				model_flux = u.Quantity(model_flux, u.erg / (u.s * u.cm**2 * u.cm))
				# BP Reading in PHOENIX model fluxes in erg/s/cm^2/cm.
				
			model_wave_file = '/d/tel2/brock/Catalogs/PHOENIX/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
			# BP File location containing wavelengths at which fluxes in PHOENIX models are reported. Not an evenly spaced grid.

			model_wavelength, wavelength_header = fits.getdata(model_wave_file, header=True)
			model_wavelength = u.Quantity(model_wavelength, u.AA)
			# BP Reading in PHOENIX model wavelengths in Angstroms.
			
			model_wavelength = self.to_air(model_wavelength)
			index = model_wavelength.argsort()
			model_wavelength, model_flux = model_wavelength[index], model_flux[index]
			# BP Converting to air wavelengths and verifying wavelengths are still in increasing order.
			model_flux = model_flux.to(u.erg / (u.s * u.cm**2 * u.AA))
			# BP Converting to Angstrom flux units.

			cs = CubicSpline(model_wavelength, model_flux)
			cs_model_flux = cs(wav_interpolate)
			# BP Creating cubic spline of spectrum model to interpolate at evenly spaced wavelengths.
			sm_model_flux = gaussian_filter(cs_model_flux, smooth_factor)
			# BP Smoothing model spectrum.

			cs_smooth = CubicSpline(wav_interpolate, sm_model_flux)
			# BP Creating cubic spline of interpolated smoothed spectrum to save as pickle. Can be called and evaluated at any wavelengths.

			### with open(pkl_file, 'wb') as f:
			### 	pickle.dump(cs_smooth, f)
				# BP Saving smoothed spectrum to pickle.

		else:
			with open(pkl_file, 'rb') as f:
				cs_smooth = pickle.load(f)
				# BP Loading the smoothed cubic spline spectrum if it exists.

		return cs_smooth

	def get_tlusty_sed(self, temperature, gravity, metallicity, smooth_factor):
		"""
		Function to return smooted cubic spline TLUSTY model spectrum to evaluate at given wavelengths.

		Parameters
		----------
		temperature : integer
			Stellar effective temperature, from ~27500 to 55000 K in steps of 2500 K.
		gravity : float
			Stellar surface gravity, from 3.00 to 4.75 in steps of 0.25.
		metallicity : float
			Stellar metallicity, choose from 0.3, 0.0, -0.3, -0.7, ..., dex.
		smooth_factor : integer
			Number of pixels of gaussian smoothing. Typically 100 for PHOENIX models.

		Returns
		-------
		cs_smooth : scipy.interpolate.CubicSpline
			Smoothed TLUSTY stellar spectrum evaluatable as a cubic spline at any wavelength.

		"""
		if np.isclose(metallicity, +0.3):
			metal_str = 'C'
		elif np.isclose(metallicity, -0.0):
			metal_str = 'G'
		elif np.islcose(metallicity, -0.3):
			metal_str = 'L'
		elif np.isclose(metallicity, -0.7):
			metal_str = 'S'
		# BP Defining TLUSTY string used to describe metalliticy. In units of 2, 1, 0.5, 0.2 times solar metallicity.
			
		pkl_file = '/d/tel2/brock/Catalogs/TLUSTY/{:.0f}t{:.0f}g{:.0f}m{:.0f}.pkl'.format(temperature, gravity*100, metallicity*100, smooth_factor)
		# BP Creating string where saved pickle file is located.
		
		if not os.path.isfile(pkl_file):
			wav_start, wav_end, wav_sep = 1400, 49000, 0.1
			wav_interpolate = np.arange(wav_start, wav_end, wav_sep)	
			# BP Creating grid of wavelengths to interpolate at.
					
			model_flux_file = '/d/tel2/brock/Catalogs/TLUSTY/{}{:.0f}g{:.0f}v10.flux'.format(metal_str, temperature, gravity*100)
			model = np.loadtxt(model_flux_file)
			# BP Defining and loading TLUSTY file with model SED.
			
			model_freq = u.Quantity(model[:,0], u.Hz)
			model_flux = u.Quantity(model[:,1], u.erg / (u.s * u.cm**2 * u.Hz))
			# BP Defining correct units for imported data.
			
			model_wavelength = model_freq.to(u.AA, equivalencies = u.spectral())
			model_flux = model_flux.to(u.erg / (u.s * u.cm**2 * u.AA), equivalencies = u.spectral_density(model_wavelength))
			# BP COnverting model fluxes into correct units.
			
			ind = (model_wavelength > wav_start * u.AA) & (model_wavelength < wav_end * u.AA)
			model_wavelength, model_flux = model_wavelength[ind], model_flux[ind]
			# BP Chopping models to more managible wavelength range.
			
			model_wavelength = self.to_air(model_wavelength)
			index = model_wavelength.argsort()
			model_wavelength, model_flux = model_wavelength[index], model_flux[index]
			# BP Converting into air wavelengths and sorting to ensure they are in increaing order.

			model_wavelength, model_ind = np.unique(model_wavelength, return_index=True)
			# BP Retrieving only the unique wavelength elements, for some reason.
			
			model_wavelength = u.Quantity(model_wavelength, u.AA)
			model_flux = model_flux[model_ind]
			# BP Converting back into Angstroms and getting the same indicies in flux.

			cs = CubicSpline(model_wavelength, model_flux)
			cs_model_flux = cs(wav_interpolate)
			# BP Creating cubic spline to interpolate at unevenly spaced wavelengths.

			sm_model_flux = gaussian_filter(cs_model_flux, smooth_factor)
			# BP Smoothing interpolated cubic spline.

			cs_smooth = CubicSpline(wav_interpolate, sm_model_flux)
			# BP Creating cubic spline of smoothed model spectrum.

			with open(pkl_file, 'wb') as f:
				pickle.dump(cs_smooth, f)
				# BP Saving smoothed spectrum to pickle.

		else:
			with open(pkl_file, 'rb') as f:
				cs_smooth = pickle.load(f)
				# BP Loading the smoothed cubic spline spectrum if it exists.

		return cs_smooth
	
	def fit_ext_lsq(self, wavelength, flux, flux_err, model_flux, ext_model):
		"""
		Function to find best fit reddening model for a given observed flux and model. Does not fit for any stellar parameters, model must be calculated externally.

		Parameters
		----------
		wavelength : array_like
			Input spectrum wavelengths.
		flux : array_like
			Input spectrum flux.
		flux_err : array_like
			Input spectrum flux error.
		model_flux : array_like
			Input unreddened model flux.
		ext_model : dust_extinciton.parameter_averages.class
			Extinction model to calculate reddening parameters with. Must be a model from the dust_extinction package.

		Returns
		-------
		model_av : float
			Best fit model Av.
		model_rv : float
			Best fit model Rv.
		model_av_err : float
			Best fit Av one sigma error.
		model_rv_err : float
			Best fit Rv one sigma error.
		final_model_spectrum : array_like
			Best fit reddened model spectrum.
		model_spectrum_lower : array_like
			Best fit one sigma lower limit model spectrum.
		model_spectrum_upper : array_like
			Best fit one sigma upper limit model spectrum.

		"""
		wavelength = u.Quantity(wavelength, u.AA)
		# BP Ensuring wavelengths are in Angstroms.
		# norm = np.argwhere((wavelength > 15000*u.AA) & (wavelength < 50000*u.AA))
		# BP Setting normailzation region to the whole wavelength.

		rv_range = ext_model.Rv_range
		# BP Extracting valid Rv range from reddening law.

		# BP Creating a function in the correct format for scipy.curve_fit.
		def reddening_model(wavelength, av, rv):
			wavelength = u.Quantity(wavelength, u.AA)
			# BP Ensuring wavelengths are in Angstroms.
			
			ext_model.Rv_range = [0.1, 10.0]
			# BP Artificially inflating the Rv range to include values less than 2.
			# BP Note that all values less than Rv~2.5 are extrapolations. None of the data stars go below 2.5, except for G23, which goes to 2.3.
			
			ext = ext_model(rv)
			ext_flux = ext.extinguish(wavelength, Av = av)
			# BP Defining reddening model/law.
				
			red_flux = model_flux * ext_flux
			# BP Creating reddened stellar spectrum.

			scale_factor = np.mean(flux)/np.mean(red_flux)
			reddened_spectrum = red_flux * scale_factor
			# BP Scaling reddened spectrum to have the same total flux as observed spectrum.

			return reddened_spectrum

		popt, pcov = curve_fit(reddening_model, wavelength, flux, bounds=((0, rv_range[0]), (np.inf, rv_range[1])), sigma=flux_err, absolute_sigma=True)
		# BP Fitting reddened spectrum with defined reddening law.	
		
		model_av, model_rv = popt[0], popt[1]
		model_err = np.sqrt(np.diag(pcov))
		model_av_err, model_rv_err = model_err[0], model_err[1]
		# BP Extracting model fit parameters.
				
		final_model_spectrum = reddening_model(wavelength, model_av, model_rv)
		# BP Creating final reddened model spectrum flux.
		
		model_rv_lower = model_rv - model_rv_err
		model_rv_upper = model_rv + model_rv_err
		# BP Calculating upper and lower Rv limits from errors.
		
		if model_rv_lower < ext_model.Rv_range[0]:
			model_rv_lower = ext_model.Rv_range[0]
			
		if model_rv_upper > ext_model.Rv_range[1]:
			model_rv_upper = ext_model.Rv_range[1]
		# BP Resetting upper and lower Rv limits to be within the range of the model if they fall outside.
		
		model_spectrum_lower = reddening_model(wavelength, model_av - model_av_err, model_rv_upper)
		model_spectrum_upper = reddening_model(wavelength, model_av + model_av_err, model_rv_lower)
		# BP Using upper and lower Rv and Av estimates from best fit parameters to construct upper and lower spectrum bounds for plotting.
		
		return model_av, model_rv, model_av_err, model_rv_err, final_model_spectrum, model_spectrum_lower, model_spectrum_upper
	
	def fit_stellar_lsq(self, wavelength, flux, flux_err, smooth_factor, temp_in = None, grav_in = None, metal_in = None):
		"""
		Function to find best fit stellar parameters for an input stellar spectrum. 
		Specify any of the three stellar parameters to constrain them to that value. 
		Uses PHOENIX stellar models from the extinction.get_phoenix_sed function.

		Parameters
		----------
		wavelength : array_like
			Input wavelength array in Angstroms.
		flux : array_like
			Input flux array.
		flux_err : array_like
			Input flux error array.
		temp_in : float, optional
			Input fixed temperature. If set temperature will not be fit. The default is None.
		grav_in : float, optional
			Input fixed surface gravity. If set gravity will not be fit. The default is None.
		metal_in : float, optional
			Input fixed metallicity. If set metallicity will not be fit. The default is None.

		Returns
		-------
		model_temperature : float
			Best fit stellar effective temperature.
		model_temp_err : float
			Standard error on the best fit temperature.
		model_gravity : float
			Best fit stellar surface gravity.
		model_grav_err : float
			Standard error on the best fit gravity.
		model_metallicity : float
			Best fit stellar metallicity.
		model_metal_err : float
			Standard error on the best fit metallicity.
		final_model_spectrum : array_like
			Output best fit model spectrum sampled at the same points as the input wavelength.
		model_spectrum_lower : array_like
			The one sigma lower best fit spectrum.
		model_spectrum_upper : array_like
			The one sigma upper best fit spectrum.

		"""
		wavelength = u.Quantity(wavelength, u.AA)
		# BP Ensuring wavelengths are in Angstroms.
		# norm = np.argwhere((wavelength > 15000*u.AA) & (wavelength < 50000*u.AA))
		# BP Setting normailzation region to the whole wavelength.

		# BP Creating a function in the correct format for scipy.curve_fit.
		def stellar_model(wavelength, temperature, gravity, metallicity):
			wavelength = u.Quantity(wavelength, u.AA)
			# BP Ensuring wavelengths are in Angstroms.
								
			cs_smooth = self.get_phoenix_sed(temperature, gravity, metallicity, smooth_factor)
			model_flux = cs_smooth(wavelength)
			# BP Getting PHOENIX stellar model flux.
			
			cont_model_flux = self.continuum_normalize(wavelength, u.Quantity(model_flux, u.erg / (u.s * u.cm**2 * u.AA)), np.repeat(model_flux.mean()*0.025, len(model_flux)))
			# BP Continuum normalizing model_flux.
		
			scale_factor = np.mean(flux)/np.mean(cont_model_flux)			
			final_flux = cont_model_flux * scale_factor
			# BP Scaling reddened spectrum to have the same total flux as observed spectrum.

			return final_flux
		
		p0 = [6000, 4.0, -0.25]
		bounds=((5500., 3.5, -2.0), (7000., 4.5, +0.5))
		# BP Specifying initial guess and bounds.

		# BP Checking which parameters were not given to fit. Holding any given parameters fixed.
		if temp_in is not None and grav_in is None and metal_in is None:
			popt, pcov = curve_fit(lambda wavelength, gravity, metallicity: stellar_model(wavelength, temp_in, gravity, metallicity), wavelength, flux, sigma=flux_err, absolute_sigma=True, p0=[p0[1], p0[2]], bounds=((bounds[0][1], bounds[0][2]),(bounds[1][1], bounds[1][2])))
			# BP Fitting stellar model with non-fixed parameters.
			model_temperature, model_gravity, model_metallicity = temp_in, popt[0], popt[1]
			model_err = np.sqrt(np.diag(pcov))
			model_temp_err, model_grav_err, model_metal_err = 0, model_err[0], model_err[1]
			# BP Extracting model fit parameters.
		
		elif temp_in is None and grav_in is not None and metal_in is None:
			popt, pcov = curve_fit(lambda wavelength, temperature, metallicity: stellar_model(wavelength, temperature, grav_in, metallicity), wavelength.value, flux, sigma=flux_err, absolute_sigma=True, p0=[p0[0], p0[2]], bounds=((bounds[0][0], bounds[0][2]),(bounds[1][0], bounds[1][2])))
			# BP Fitting stellar model with non-fixed parameters.
			model_temperature, model_gravity, model_metallicity = popt[0], grav_in, popt[1]
			model_err = np.sqrt(np.diag(pcov))
			model_temp_err, model_grav_err, model_metal_err = model_err[0], 0, model_err[1]
			# BP Extracting model fit parameters.
		
		elif temp_in is None and grav_in is None and metal_in is not None:
			popt, pcov = curve_fit(lambda wavelength, temperature, gravity: stellar_model(wavelength, temperature, gravity, metal_in), wavelength.value, flux, sigma=flux_err, absolute_sigma=True, p0=[p0[0], p0[1]], bounds=((bounds[0][0], bounds[0][1]),(bounds[1][0], bounds[1][1])))
			# BP Fitting stellar model with non-fixed parameters.
			model_temperature, model_gravity, model_metallicity = popt[0], popt[1], metal_in
			model_err = np.sqrt(np.diag(pcov))
			model_temp_err, model_grav_err, model_metal_err = model_err[0], model_err[1], 0
			# BP Extracting model fit parameters.
		
		elif temp_in is not None and grav_in is not None and metal_in is None:
			popt, pcov = curve_fit(lambda wavelength, metallicity: stellar_model(wavelength, temp_in, grav_in, metallicity), wavelength.value, flux, sigma=flux_err, absolute_sigma=True, p0=[p0[2]], bounds=((bounds[0][2]),(bounds[1][2])))
			# BP Fitting stellar model with non-fixed parameters.
			model_temperature, model_gravity, model_metallicity = temp_in, grav_in, popt[0]
			model_err = np.sqrt(np.diag(pcov))
			model_temp_err, model_grav_err, model_metal_err = 0, 0, model_err[0]
			# BP Extracting model fit parameters.
		
		elif temp_in is not None and grav_in is None and metal_in is not None:
			popt, pcov = curve_fit(lambda wavelength, gravity: stellar_model(wavelength, temp_in, gravity, metal_in), wavelength.value, flux, sigma=flux_err, absolute_sigma=True, p0=[p0[1]], bounds=((bounds[0][1]),(bounds[1][1])))
			# BP Fitting stellar model with non-fixed parameters.
			model_temperature, model_gravity, model_metallicity = temp_in, popt[0], metal_in
			model_err = np.sqrt(np.diag(pcov))
			model_temp_err, model_grav_err, model_metal_err = 0, model_err[0], 0
			# BP Extracting model fit parameters.
		
		elif temp_in is None and grav_in is not None and metal_in is not None:
			popt, pcov = curve_fit(lambda wavelength, temperature: stellar_model(wavelength, temperature, grav_in, metal_in), wavelength.value, flux, sigma=flux_err, absolute_sigma=True, p0=[p0[0]], bounds=((bounds[0][0]),(bounds[1][0])))
			# BP Fitting stellar model with non-fixed parameters.
			model_temperature, model_gravity, model_metallicity = popt[0], grav_in, metal_in
			model_err = np.sqrt(np.diag(pcov))
			model_temp_err, model_grav_err, model_metal_err = model_err[0], 0, 0
			# BP Extracting model fit parameters.
		
		elif temp_in is not None and grav_in is not None and metal_in is not None:
			model_temperature, model_gravity, model_metallicity = temp_in, grav_in, metal_in
			model_temp_err, model_grav_err, model_metal_err = 0, 0, 0
			# BP Extracting model fit parameters.
		
		else:
			popt, pcov = curve_fit(stellar_model, wavelength.value, flux, sigma=flux_err, absolute_sigma=True, p0=p0, bounds=bounds)
			# BP Fitting reddened spectrum with defined reddening law.	
			model_temperature, model_gravity, model_metallicity = popt[0], popt[1], popt[2]
			model_err = np.sqrt(np.diag(pcov))
			model_temp_err, model_grav_err, model_metal_err = model_err[0], model_err[1], model_err[2]
			# BP Extracting model fit parameters.
							
		final_model_spectrum = stellar_model(wavelength, model_temperature, model_gravity, model_metallicity)
		# BP Creating final reddened model spectrum flux.
		
		#metal_up = 
		
		model_spectrum_lower = 0#stellar_model(wavelength, model_temperature - model_temp_err, model_gravity + model_grav_err, model_metallicity + model_metal_err)
		model_spectrum_upper = 0#stellar_model(wavelength, model_temperature + model_temp_err, model_gravity - model_grav_err, model_metallicity - model_metal_err)
		# BP Using upper and lower Rv and Av estimates from best fit parameters to construct upper and lower spectrum bounds for plotting.
		# UPPER AND LOWER MIGHT BE WRONG
		
		return model_temperature, model_temp_err, model_gravity, model_grav_err, model_metallicity, model_metal_err, final_model_spectrum, model_spectrum_lower, model_spectrum_upper
	
	def to_air(self, wave_vac):
		"""
		Function to convert from vacuum wavelengths to air wavelengths at STP to match HST observations.

		Parameters
		----------
		wave_vac : array
			Input array of vacuum wavelengths in Angstroms.

		Returns
		-------
		wave_air : array
			Output array of air wavelengths in Angstroms.
			
		"""
		wave_vac = u.Quantity(wave_vac, u.AA)
		# BP Ensuring input wavelengths are in Angstroms.
		s = 1/wave_vac.to(u.um)
		# BP Converting to inverse microns. 
		n_stp = 1.0 + 0.0000834254 + 0.02406147*(u.um**-2) / (130*(u.um**-2) - s**2) + 0.00015998*(u.um**-2) / (38.9*(u.um**-2) - s**2)
		# BP Calculating index of refraction at STP as per http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
		wave_air = wave_vac/n_stp
		# BP Converting into air wavelengths.
		
		return wave_air
	
	def plot_best_fit(self, name, wavelength, flux, flux_err, model_flux, model_flux_lower, model_flux_upper, redchisq, temperature, gravity, metallicity, ext_model, av, rv, av_err, rv_err):
		"""
		Comprehensive function to plot raw spectra along with the best fit reddening model. All calculations occur externally, input final analysis.

		Parameters
		----------
		name : string
			Star name to plot.
		wavelength : array_like
			Input wavelength array in Angstroms.
		flux : array_like
			Input flux array.
		flux_err : array_lik
			Input flux error array.
		model_flux : array_like
			Input model flux.
		model_flux_lower : array_like
			Input lower limit for model flux using rv and av errors.
		model_flux_upper : array_like
			Input upper limit for model flux using rv and av errors.
		redchisq : integer
			Reduced chi squared of model.
		temperature : string
			Stellar effective temperature.
		metallicity : string
			Stellar metallicity.
		gravity : string
			Stellar surface gravity.
		ext_model : string
			Extinction model used to calculate reddening parameters.
		av : float
			Calculated model Av.
		rv : float
			Calculated model Rv.
		av_err : float
			Model Av error.
		rv_err : float
			Model Rv error.

		Returns
		-------
		None.

		"""
		fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(8,6), gridspec_kw={'height_ratios': [2.5, 1], 'wspace':0,'hspace':0}, sharex='all', layout='tight')
		# BP Initializing plotting axes for two plot figure with residuals.

		diff_lower = np.log10(flux - flux_err)
		diff_upper = np.log10(flux + flux_err)
		diff_lower[np.isnan(diff_lower)] = 0
		# BP Calculating errors in log space.

		yerr_lower = abs(np.log10(flux) - diff_lower)
		yerr_upper = abs(np.log10(flux) - diff_upper)
		# BP Translating log errors into upper and lower error bars.

		ax1.plot(np.log10(wavelength), np.log10(model_flux), color = 'r', label = 'Best Fit Reddened Model', zorder=5)
		ax1.plot(np.log10(wavelength), np.log10(flux), color='b', linestyle='-')
		ax1.errorbar(np.log10(wavelength), np.log10(flux), yerr=[yerr_lower, yerr_upper], color='b', linestyle='',label='HST Spectra + Photometry', capsize=2, zorder=10)
		ax1.fill_between(np.log10(wavelength), np.log10(model_flux_lower), np.log10(model_flux_upper), alpha = 0.25, facecolor='r')
		# BP Plotting model spectrum with error regions and data.

		ax2.set_xlabel('log Wavelength (Ang.)')
		ax1.tick_params(axis='both', direction='in', which='both')
		ax1.legend(loc='lower right')
		ax1.set_ylabel('log F$_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ Ang$^{-1}$) ')
		# BP Setting labels.

		residuals = flux - model_flux
		sigmas = residuals/flux_err
		upper,lower = sigmas.max(),sigmas.min()
		# BP Calculating residuals in terms of sigmas and getting max and min for plotting.
		
		ax2.set_ylabel('Residuals ($\sigma$)')
		ax2.plot(np.log10(wavelength), sigmas, color = 'k', linestyle='', marker = '.', markersize=5, label='$\chi^2_{{red}}:$ {:.3f}'.format(redchisq))
		ax2.axhline(0,color = 'k', linestyle='--', alpha = 0.25,linewidth=1, dashes = (5,5))
		ax2.set_ylim([lower*1.25,upper*1.25])
		ax2.legend(loc = 'upper right', handlelength=0, handletextpad=0, fancybox=True, markerscale=1e-20)
		ax2.tick_params(axis='both', direction='in', which='both', top=True)
		# BP Setting plotting axes and parameters.
				
		ax1.annotate('{}:\n$T_{{eff}}$={:.1f}, [Fe/H]={:.3f}, log $g$={:.3f}\n{}:\nA$_{{V}}$={:.3f}±{:.3f}, R$_{{V}}$={:.3f}±{:.3f}'.format(name.replace('pk','PK-'), temperature, metallicity, gravity, ext_model, av, av_err, rv, rv_err), size=18, xy=(0.965,0.535), xycoords='figure fraction', xytext=(-20,-10), textcoords='offset points', ha='right', va='bottom')
		plt.savefig("/d/users/brock/paper_figures/{}_{}.pdf".format(name, ext_model))
		plt.show()
		# BP Plotting text and saving plot for PK targets.

# =============================================================================
# 		ax1.annotate('{}:\n$T_{{eff}}$={}, [Fe/H]=${}$, log $g$={}\n{}:\nA$_{{V}}$={:.3f}±{:.3f}, R$_{{V}}$={:.3f}±{:.3f}'.format(name, temperature, metallicity, gravity, ext_model, av, av_err, rv, rv_err), size=18, xy=(0.21, 0.4), xycoords='figure fraction', xytext=(-20,-10), textcoords='offset points', ha='left', va='bottom')
# 		plt.savefig("/d/users/brock/paper_figures/{}_{}.pdf".format(name, ext_model))
# 		plt.show()
# 		# BP Plotting text and saving plot for Zeta Oph.
# =============================================================================
		
	def chi_squared(self, data, model, err, deg=3):
		"""
		Function to calculate chi squared and reduced chi squared of a model and data.

		Parameters
		----------
		data : array_like
			Input observed data.
		model : array_like
			Input model data.
		err : array_like
			Input observed data error.
		deg : integer, optional
			The degrees of freedom of the model. The default is 3 for Av, Rv, and normalization.

		Returns
		-------
		chisq : integer
			Chi squared of model.
		red_chisq : integer
			Reduced chi squared of data.

		"""
		chisq = np.sum(((data-model)/err)**2)
		# BP Calculating chi squared of the model.
		n = data.size - 1 - deg
		# BP Calculating degrees for the model.
		
		red_chisq = chisq/n
		# BP Calculating reduced chi squared.

		return chisq, red_chisq
	
	def continuum_normalize(self, wavelength, flux, flux_err, plot=False):
		"""
		Function to continuum normalize a given input spectrum using specutils. Normalizes to 1.

		Parameters
		----------
		wavelength : array_like
			Input spectrum wavelengths.
		flux : array_like
			Input spectrum flux.
		flux_err : array_like
			Input flux error array.
		plot : boolean, optional
			Flag for plotting of continuum fit. Default is False.

		Returns
		-------
		cont_norm_spec : array_like
			Continuum normalized spectrum flux. No units.

		"""	
		wavelength = u.Quantity(wavelength, u.AA).value
		flux = u.Quantity(flux, u.erg / (u.s * u.cm**2 * u.AA)).value
		flux_err = u.Quantity(flux_err, u.erg / (u.s * u.cm**2 * u.AA)).value
		# BP Ensuring all inputs are in the correct units.
		
		ind = (wavelength == wavelength)
		#ind = find_peaks(flux)[0]
		# BP Finding just the peaks (local maxima) to fit with the continuum.
		wavelength_fit, flux_fit, flux_err_fit = wavelength[ind], flux[ind], flux_err[ind]
		# BP Indexing at only the peaks.
		
		cheby = models.Chebyshev1D(5)
		fit = fitting.LMLSQFitter()
		# BP Defining the model to fit and the fitting model.
		
		or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip, niter=10, sigma_lower = 2, sigma_upper = 5.0)
		# Fitting the continuum with a low order chebyshev polynomial using iterative fitting with outlier removal.
		continuum, mask = or_fit(cheby, wavelength_fit, flux_fit, weights=1.0/flux_err_fit)
		filter_continuum = np.ma.masked_array(flux_fit, mask=mask)
		# BP Evaluating the fit at given wavelengths and finding outliers.
		
		cont_norm_flux = flux / (continuum(wavelength))
		# BP Calculating the flux after continuum fitting.
	
		if plot==True:
			fig, ax = plt.subplots(figsize=(10,4), layout='tight')
			# BP Creating plotting axes.
		
			plt.plot(wavelength, continuum(wavelength), linestyle = '-', marker = '', color = 'k')
			plt.plot(wavelength[ind], flux_fit, marker = '.', linestyle = '', color = 'r')
			plt.plot(wavelength[ind], filter_continuum, marker = 'o', linestyle = '', color = 'g')
			# BP Plotting best fit continuum and point rejection.
	
			plt.show()

		return cont_norm_flux

	def interpolate_phoenix(self):
		"""
		Function to interpolate PHOENIX stellar models using Trilinear interpolation. 

		Parameters
		----------
		None.

		Returns
		-------
		interp : scipy.RegularGridInterpolator
			Interpolated PHOENIX grid spectra that can be individually called at specific stellar parameters.

		"""
		interp_file = '/d/tel2/brock/Catalogs/PHOENIX/interp_grid.pkl'
		# BP Defining location of pickle file continaing interpolating polynomial.
		
		if not os.path.isfile(interp_file):
			temps = np.arange(5000, 6600, 100)
			# BP Creating temperature grid covering suspected temperatures of targets.
			
			metals = np.array([-1.0, -0.5, -0.0, +0.5])
			# BP Creating metallicity grid at -2.0, -1.0, -0.5, -0.0, and +0.5 dex.
			
			gravs = np.array([4.0, 4.5])
			# BP Creating gravity grid from 1.0 to 6.0 every 0.5 dex.
					
			flux_arr = np.zeros((len(temps), len(gravs), len(metals), 1569128))
			# BP Creating zeros array to write all spectra into to input into interpolation function.
			
			for i in range(len(temps)):
				for j in range(len(gravs)):
					for k in range(len(metals)):
						temp = temps[i]
						grav = gravs[j]
						metal = metals[k]
						# BP Looping over all temperatures, metallicities, and gravities defined above.
						
						flux_file = '/d/tel2/brock/Catalogs/PHOENIX/lte{:05d}-{:.2f}{:+.1f}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'.format(temp, grav, metal)
						flux, header = fits.getdata(flux_file, header=True)
						flux_arr[i,j,k, :] = flux
						# BP Extracting the flux from each of the ~90 corresponding models and saving to one large array.
			
			interp = RegularGridInterpolator((temps, gravs, metals), flux_arr)
			# BP Writing out interpolation as callable function.
			
			with open(interp_file, 'wb') as f:
				pickle.dump(interp, f)
				# BP Saving interpolation function as pickle file.
		
		else:
			with open(interp_file, 'rb') as f:
				interp = pickle.load(f)
				# BP Loading interpolation 
		
		return interp
		