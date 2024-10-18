'''Importable script to calculate reddening parameters from given stellar spectra. Call as import reddening.'''

import numpy as np
import pickle
import os
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
from scipy.interpolate import CubicSpline, interp1d
from scipy.ndimage import gaussian_filter
from scipy import signal
from scipy.optimize import curve_fit

from dust_extinction.parameter_averages import F19, CCM89, VCG04, G23
# BP Neccessary imports.

class spectrum(object):
# =============================================================================
# 	def to_air(self, wave_vac):
# 		'''
# 		Converts input wavelengths from space wavelengths to air wavelengths.
# 		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 		Inputs:
# 			wave_vac: astropy.units.quantity.Quantity. Array of input vacuum wavelengths.
# 		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 		Outputs:
# 			wave_air: astropy.units.quantity.Quantity. Array of converted air wavelengths.
# 		'''
# 		wave_vac = u.Quantity(wave_vac, u.AA)
# 		# BP Ensuring input wavelengths are in Angstroms.
# 		s = 1/wave_vac.to(u.um)
# 		# BP Converting to inverse microns. 
# 		n_stp = 1.0 + 0.0000834254 + 0.02406147*(u.um**-2) / (130*(u.um**-2) - s**2) + 0.00015998*(u.um**-2) / (38.9*(u.um**-2) - s**2)
# 		# BP Calculating index of refraction at STP as per http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
# 		wave_air = wave_vac/n_stp
# 		# BP Converting into air wavelengths.
# 		
# 		return wave_air
# =============================================================================

# =============================================================================
# 	def __init__(self, wavelengths, flux, flux_err, temperature, gravity, metallicity, extinction_1, reddening_1, extinction_2, reddening_2, star, ext_model, smooth_factor):
# 		'''Initializes all global variables that are used throughout the script.
# 		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 		Inputs:
# 			wavelengths: astropy.units.quantity.Quantity. Array of observed wavelengths.
# 			flux: astropy.units.quantity.Quantity. Array of observed fluxes.
# 			flux_error: astropy.units.quantity.Quantity. Array of standard deviations for flux measurements.
# 			temperature: 
# 			
# 		'''
# 		#### Add np.sort to automatically place wavelengths (attached with fluxes+errors) in increasing order
# 		self._wavelengths = wavelengths
# 		# self._wavelengths = self.to_air(wavelengths)
# 		self._flux = flux
# 		self._flux_err = flux_err
# 		self._temperature = temperature
# 		self._gravity = gravity
# 		self._metallicity = metallicity
# 		self._extinction_1 = extinction_1
# 		self._reddening_1 = reddening_1
# 		self._extinction_2 = extinction_2
# 		self._reddening_2 = reddening_2
# 		self._star = star
# 		self._ext_model = ext_model
# 		self._smooth_factor = smooth_factor
# =============================================================================


# =============================================================================
# 	def CCM89(self, microns, Av, Rv):
# 		inverse = np.array(1./microns)
# 		IR = np.array(np.nonzero(inverse < 1.1)).ravel()
# 		Optical = np.array(np.nonzero((inverse >= 1.1) & (inverse <= 3.3))).ravel()
# 		UV = np.array(np.nonzero((inverse > 3.3) & (inverse <= 8.))).ravel()
# 
# 		if len(IR) > 0:
# 			a_IR = 0.574 * (inverse[IR]**1.61)
# 			b_IR = -0.527 * (inverse[IR]**1.61)
# 			lam_IR = (a_IR + b_IR/Rv) * Av
# 		else:
# 			lam_IR = np.array([])
# 
# 		if len(Optical) > 0:
# 			y = inverse[Optical]-1.82
# 			a_opt = 1. + 0.17699*y - 0.50447*(y**2) - 0.02427*(y**3) + 0.72085*(y**4) + 0.01979*(y**5) - 0.77530*(y**6) + 0.32999*(y**7)
# 			b_opt = 1.41338*y + 2.28305*(y**2) + 1.07233*(y**3) - 5.38434*(y**4) - 0.62251*(y**5) + 5.30260*(y**6) - 2.09002*(y**7)
# 			lam_opt=(a_opt + b_opt/Rv) * Av
# 		else:
# 			lam_opt=np.array([])
# 		
# 		if len(UV) > 0:
# 			UV_2 = np.array(np.nonzero((inverse >= 5.9) & (inverse <= 8.))).ravel()
# 			UV_1 = np.array(np.nonzero((inverse < 5.9) & (inverse > 3.3))).ravel()
# 
# 			Fa_2 = -0.04473*(inverse[UV_2]-5.9)**2 - 0.009779*(inverse[UV_2]-5.9)**3
# 			Fb_2 = 0.2130*(inverse[UV_2]-5.9)**2 + 0.1207*(inverse[UV_2]-5.9)**3
# 			
# 			Fa_1 = UV_1 - UV_1	
# 			Fb_1 = Fa_1		
# 
# 			F_a = np.concatenate([Fa_2, Fa_1])
# 			F_b = np.concatenate([Fb_2, Fb_1])
# 
# 			a_UV = 1.752 - 0.316*inverse[UV] - 0.104/((inverse[UV]-4.67)**2 + 0.341) + F_a
# 			b_UV = -3.090 + 1.825*inverse[UV] + 1.206/((inverse[UV]-4.62)**2 + 0.263) + F_b
# 	   
# 			lam_UV = (a_UV + b_UV/Rv) * Av
# 		else:
# 			lam_UV=np.array([])
# 
# 		A_lam = np.concatenate([lam_UV,lam_opt,lam_IR])
# 
# 		return A_lam
# =============================================================================


# =============================================================================
# 	def VCG04(self, microns, Av, Rv):
# 		inverse = np.array(1./microns)
# 		IR = np.array(np.nonzero(inverse < 1.1)).ravel()
# 		Optical = np.array(np.nonzero((inverse >= 1.1) & (inverse <= 3.3))).ravel()
# 		UV = np.array(np.nonzero((inverse > 3.3) & (inverse <= 8.))).ravel()
# 
# 		if len(IR) > 0:
# 			a_IR = 0.574 * (inverse[IR]**1.61)
# 			b_IR = -0.527 * (inverse[IR]**1.61)
# 			lam_IR = (a_IR + b_IR/Rv) * Av
# 		else:
# 			lam_IR = np.array([])
# 
# 		if len(Optical) > 0:
# 			y = inverse[Optical]-1.82
# 			a_opt = 1. + 0.17699*y - 0.50447*(y**2) - 0.02427*(y**3) + 0.72085*(y**4) + 0.01979*(y**5) - 0.77530*(y**6) + 0.32999*(y**7)
# 			b_opt = 1.41338*y + 2.28305*(y**2) + 1.07233*(y**3) - 5.38434*(y**4) - 0.62251*(y**5) + 5.30260*(y**6) - 2.09002*(y**7)
# 			lam_opt=(a_opt + b_opt/Rv) * Av
# 		else:
# 			lam_opt=np.array([])
# 		
# 		if len(UV) > 0:
# 			UV_2 = np.array(np.nonzero((inverse >= 5.9) & (inverse <= 8.))).ravel()
# 			UV_1 = np.array(np.nonzero((inverse < 5.9) & (inverse > 3.3))).ravel()
# 
# 			Fa_2 = -0.0077*(inverse[UV_2]-5.9)**2 - 0.0030*(inverse[UV_2]-5.9)**3
# 			Fb_2 = 0.2060*(inverse[UV_2]-5.9)**2 + 0.0550*(inverse[UV_2]-5.9)**3
# 			
# 			Fa_1 = UV_1 - UV_1	
# 			Fb_1 = Fa_1		
# 
# 			F_a = np.concatenate([Fa_2, Fa_1])
# 			F_b = np.concatenate([Fb_2, Fb_1])
# 
# 			#a_UV = 1.708 - 0.215*inverse[UV] - 0.134/((inverse[UV]-4.558)**2 + 0.566) + F_a		
# 			# changing a_UV original:
# 			a_UV = 1.808 - 0.215*inverse[UV] - 0.134/((inverse[UV]-4.558)**2 + 0.566) + F_a
# 			b_UV = -2.350 + 1.403*inverse[UV] + 1.103/((inverse[UV]-4.587)**2 + 0.263) + F_b
# 	   
# 			lam_UV = (a_UV + b_UV/Rv) * Av
# 		else:
# 			lam_UV=np.array([])
# 
# 		A_lam = np.concatenate([lam_UV,lam_opt,lam_IR])
# 
# 		return A_lam
# =============================================================================

# =============================================================================
# 	def chisquared(self, data, model,deg=2,err=None):
# 		chisq=np.sum(((data-model)/err)**2)  
# 		n=data.size - 1 - deg
# 
# 		return chisq, chisq/n
# =============================================================================


# =============================================================================
# 	def plot_sed(self, model_spectrum, redchisq):
# 		plt.rc('axes', labelsize=16)
# 		plt.rc('figure', titlesize=16)
# 		plt.rc('xtick', labelsize=16)
# 		plt.rc('ytick', labelsize=16)
# 		plt.rc('legend', fontsize=16)
# 
# 		fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(8,6), gridspec_kw={'height_ratios': [2.5, 1], 'wspace':0,'hspace':0}, sharex='all', layout='tight')
# 
# 		#ax1.plot(self._wavelengths, model_spectrum, color = 'r', label = 'PHOENIX Model', zorder=5)
# 		#ax1.plot(self._wavelengths, self._flux, color='b', linestyle='-')
# 		#ax1.errorbar(self._wavelengths, self._flux, yerr=self._flux_err, color='b', linestyle='',label='STIS Spectra + Photometry', capsize=2)
# 		#ax1.set_yscale('log')
# 
# 		diff_lower = np.log10(self._flux - self._flux_err)
# 		diff_upper = np.log10(self._flux + self._flux_err)
# 		diff_lower[np.isnan(diff_lower)] = 0
# 
# 		yerr_lower = abs(np.log10(self._flux) - diff_lower)
# 		yerr_upper = abs(np.log10(self._flux) - diff_upper)
# 
# 		#ax1.plot(np.log10(self._wavelengths), np.log10(model_spectrum), color = 'r', label = 'MT-NextGen Model', zorder=5)
# 		ax1.plot(np.log10(self._wavelengths), np.log10(model_spectrum), color = 'r', label = 'PHOENIX Model', zorder=5)
# 		ax1.plot(np.log10(self._wavelengths), np.log10(self._flux), color='b', linestyle='-')
# 		#ax1.errorbar(np.log10(self._wavelengths), np.log10(self._flux), yerr=[yerr_lower, yerr_upper], color='b', linestyle='',label='IUE Spectra + Photometry', capsize=2)
# 		ax1.errorbar(np.log10(self._wavelengths), np.log10(self._flux), yerr=[yerr_lower, yerr_upper], color='b', linestyle='',label='HST Spectra + Photometry', capsize=2)
# 
# 
# 		ax2.set_xlabel('log Wavelength (Ang.)')
# 		ax1.tick_params(axis='both', direction='in', which='both')
# 		ax1.legend(loc='lower right')
# 		ax1.set_ylabel('log F$_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ Ang$^{-1}$) ')
# 
# 		max_spec = self._flux.max()*1.1
# 
# 		#ax1.set_xlim([1400,22000])
# 		#ax1.set_ylim(bottom = -5.75)
# 
# 		residuals = self._flux - model_spectrum
# 		sigmas = residuals/self._flux_err
# 		upper,lower = sigmas.max(),sigmas.min()
# 		
# 		ax2.set_ylabel('Residuals ($\sigma$)')
# 		ax2.plot(np.log10(self._wavelengths), sigmas, color = 'k', linestyle='', marker = '.', markersize=5, label='$\chi^2_{{red}}:$ {:.3f}'.format(redchisq))
# 		ax2.axhline(0,color = 'k', linestyle='--', alpha = 0.25,linewidth=1, dashes = (5,5))
# 		ax2.set_ylim([lower*1.25,upper*1.25])
# 		ax2.legend(loc = 'upper right', handlelength=0, handletextpad=0, fancybox=True, markerscale=1e-20)
# 		ax2.tick_params(axis='both', direction='in', which='both', top=True)
# 
# 		#ax1.text(3.3, 2.5, '$\zeta$ Ophiuchi:\n$T_{{eff}}$={}, [Fe/H]={}, log $g$={}\n{}:\nA$_{{V}}$={:.3f}, R$_{{V}}$={:.3f}'.format(self._temperature, self._metallicity, self._gravity, self._ext_model, self._extinction_1, self._reddening_1, self._extinction_2, self._reddening_2), size = 18) #, A$_{{V2}}$={:.3f}, R$_{{V2}}$={:.3f}
# 		ax1.text(4.36, -3.4, 'PK2021-27:\n$T_{{eff}}$={}, [Fe/H]={}, log $g$={}\n{}:\nA$_{{V}}$={:.3f}, R$_{{V}}$={:.3f}'.format(self._temperature, self._metallicity, self._gravity, self._ext_model, self._extinction_1, self._reddening_1, self._extinction_2, self._reddening_2), size = 18, horizontalalignment='right') #, A$_{{V2}}$={:.3f}, R$_{{V2}}$={:.3f}
# 		plt.savefig(self._star + "_SED.pdf")
# 		plt.show()
# =============================================================================

# =============================================================================
# 	def plot_raw(self, temps, limits=None):
# 		plt.rc('axes', labelsize=16)
# 		plt.rc('figure', titlesize=16)
# 		plt.rc('xtick', labelsize=16)
# 		plt.rc('ytick', labelsize=16)
# 		plt.rc('legend', fontsize=16)
# 		
# 		fig, ax = plt.subplots(figsize=(6,8), layout='tight')
# 		stars = self._star[::-1]
# 		temps = temps[::-1]
# 
# 		for i in range(len(stars)):
# 			star = stars[i]
# 			basedir = '/d/tel2/brock/Data/HST/cyc29/data/'
# 			s1w,s1f,s1e=np.loadtxt(basedir + star + '.dat',unpack=True,skiprows=1)
# 			s1e = np.array(s1e)#/1e-14
# 			s1f = np.array(s1f)#/1e-14
# 			s1f[s1f < 0] = 0
# 			
# 			s1f = gaussian_filter(s1f, 1)
# 
# 			if limits is not None:
# 				wav = limits[i]
# 				s1w, s1e, s1f = s1w[s1w>wav], s1e[s1w>wav], s1f[s1w>wav]
# 				good = np.where((s1w > 1400) & (s1w < 3155))
# 				self._wavelengths, self._flux_err, self._flux = s1w[good], s1e[good], s1f[good]
# 			else:
# 				good = np.where((s1w > 1400) & (s1w < 3155))
# 				self._wavelengths, self._flux_err, self._flux = s1w[good], s1e[good], s1f[good]
# 			
# 			shift = 1.75
# 
# 			ax.plot(self._wavelengths, np.log10(self._flux) + shift*i, linestyle='-', color='b')
# 			ax.text(1790, -1.15+shift*i, 'PK-' + star[-2:] + ', $T_{eff}$=' + temps[i], size = 18, horizontalalignment='left')
# 				
# 		ax.axvline(1931, color='k', linestyle=':', linewidth = 0.75)
# 		ax.text(1945, -2.5, 'C I', rotation=90, fontsize=14)
# 		
# 		ax.axvline(2745, color='k', linestyle=':', linewidth = 0.75)
# 		ax.text(2695, -2.5, 'Fe II/Fe I', rotation=90, fontsize=14)
# 		
# 		ax.axvline(2799, color='k', linestyle=':', linewidth = 0.75)
# 		ax.text(2755, -2.5, 'Mg II', rotation=90, fontsize = 14)
# 		
# 		ax.axvline(2852, color='k', linestyle=':', linewidth = 0.75)
# 		ax.text(2805, -2.5, 'Mg I', rotation=90, fontsize=14)
# 		
# 		ax.axvline(2881, color='k', linestyle=':', linewidth = 0.75)
# 		ax.text(2890, -2.5, 'Si I', rotation=90, fontsize=14)
# 
# 		ax.set_xlabel('Wavelength (Ang.)')
# 		ax.tick_params(axis='both', direction='in', which='both')
# 		ax.set_ylabel('log F$_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ Ang$^{-1}$) + offset')
# 		plt.show()
# 		
# =============================================================================
# =============================================================================
# 	def snr(self):
# 		plt.rc('axes', labelsize=16)
# 		plt.rc('figure', titlesize=16)
# 		plt.rc('xtick', labelsize=16)
# 		plt.rc('ytick', labelsize=16)
# 		plt.rc('legend', fontsize=16)
# 		
# 		fig, ax = plt.subplots(figsize=(8,5), layout='tight')
# 		stars = self._star[::-1]
# 
# 		snr_lims = []
# 
# 		for i in range(len(stars)):
# 			star = stars[i]
# 			basedir = '/d/tel2/brock/Data/HST/cyc29/data/'
# 			s1w,s1f,s1e=np.loadtxt(basedir + star + '.dat',unpack=True,skiprows=1)
# 			s1e = np.array(s1e)#/1e-14
# 			s1f = np.array(s1f)#/1e-14
# 			
# 			good = np.where((s1w > 1400) & (s1w < 3200))
# 			self._wavelengths, self._flux_err, self._flux = s1w[good], s1e[good], s1f[good]
# 			
# 			run_avg, run_sig = np.array([]), np.array([])
# 			
# 			width = 25
# 			
# 			for k in range(width,len(self._wavelengths)-width):
# 				run_avg = np.append(run_avg, np.mean(self._flux[k:width*2+k]))
# 				run_sig = np.append(run_sig, np.std(self._flux[k:width*2+k]))
# 				
# 			if i==0:
# 				snr = np.array([run_avg/run_sig])
# 			else:
# 				snr = np.vstack((snr, run_avg/run_sig))
# 			
# 			snr_lims.append(self._wavelengths[width:-width][np.where(snr[i] < 2)[0][-1]])
# 									
# 			ax.plot(self._wavelengths[width:-width], snr[i], linestyle='-', marker='', label=star[-2:])
# 			ax.axhline(2, linestyle='--', color='gray')
# 			
# 		ax.axvline(1931, color='k', linestyle=':', linewidth = 0.75)
# 		ax.text(1931, 15, 'C I', rotation=90)
# 		
# 		ax.axvline(2745, color='k', linestyle=':', linewidth = 0.75)
# 		ax.text(2745, 15, 'Fe II/Fe I', rotation=90)
# 		
# 		ax.axvline(2799, color='k', linestyle=':', linewidth = 0.75)
# 		ax.text(2799, 15, 'Mg II', rotation=90)
# 		
# 		ax.axvline(2852, color='k', linestyle=':', linewidth = 0.75)
# 		ax.text(2852, 15, 'Mg I', rotation=90)
# 		
# 		ax.axvline(2881, color='k', linestyle=':', linewidth = 0.75)
# 		ax.text(2881, 15, 'Si I', rotation=90)
# 			
# 		ax.tick_params(axis='both', direction='in', which='both')
# 		ax.set_ylabel('SNR')
# 		ax.set_xlabel('Wavelength (Ang.)')
# 		ax.legend(loc='upper left')
# 		plt.show()
# 			
# 		return snr_lims
# =============================================================================


# =============================================================================
# 	def fit_sed(self, plot=0):
# 		pkl_file = '/d/tel2/brock/Catalogs/PHOENIX/lte' + self._temperature.zfill(5) + "-" + self._gravity + self._metallicity + '-' + str(self._smooth_factor) + ".PHOENIX-ACES-AGSS-COND-2011-HiRes.pkl"
# 		if not os.path.isfile(pkl_file):
# 			lam1, lam2, cdelt1 = 1400, 49000, 0.1
# 
# 			lamA1 = np.arange(lam1,lam2,cdelt1)	
# 
# 			if (self._star == 'zetaOph_IUE') | (self._star == 'sedzetaOph'):
# 				wavs,flux=np.loadtxt("/d/zem1/hak/chip/Ashley/Reddening/zeta_BT-NextGenmodel.dat",unpack=True,skiprows=100,max_rows=363000)  #31000K AGSS2009 MT-NextGen model, logg=4, flux in erg/s/cm^2/Ang, wavs in Ang
# 				wavs=wavs[0:-1:7]
# 				flux=flux[0:-1:7]
# 
# 			else:
# 				wavefile="/d/tel2/brock/Catalogs/PHOENIX/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
# 
# 				fluxfile = '/d/tel2/brock/Catalogs/PHOENIX/lte' + self._temperature.zfill(5) + "-" + self._gravity + self._metallicity + ".PHOENIX-ACES-AGSS-COND-2011-HiRes_test.fits"
# 				flux, header = fits.getdata(fluxfile, header=True) 	# fluxes in erg/s/cm^2/cm
# 				wavs, header2 = fits.getdata(wavefile, header=True) 	# wavelengths in Ang
# 				wavs = self.to_air(wavs)				# convert to air wavelengths
# 				index = wavs.argsort()					# make sure wavs array is still in increasing order
# 				wavs, flux = wavs[index], flux[index]
# 				flux = flux*1e-8
# 
# 			cs = CubicSpline(wavs, flux)
# 			newflux = cs(lamA1)
# 			smflux = gaussian_filter(newflux, self._smooth_factor)
# 
# 			cs2 = CubicSpline(lamA1, smflux)
# 
# 			with open(pkl_file, 'wb') as f:
# 				pickle.dump(cs2, f)			
# 
# 			finalflux = cs2(self._wavelengths)
# 		else:
# 			with open(pkl_file, 'rb') as f:
# 				cs2 = pickle.load(f)
# 			finalflux = cs2(self._wavelengths)
# 
# 		#self._flux = gaussian_filter(self._flux, self._smooth_factor/20)
# 
# 		#norm = np.argwhere((self._wavelengths < 50000) & (self._wavelengths > 30000))
# 		norm = np.argwhere((self._wavelengths < 50000) & (self._wavelengths > 1600))
# 
# 		#if self._ext_model == 'CCM89':
# 			#Extinc1 = CCM89(self._reddening_1)
# 			#Extinc1F = Extinc1.extinguish(1/self._wavelengths*1e4, Av =  self._extinction_1)
# 			#Extinc2 = CCM89(self._reddening_2)
# 			#Extinc2F = Extinc2.extinguish(1/self._wavelengths*1e4, Av = self._extinction_2)
# 
# 			#Extflux1 = finalflux*Extinc1F
# 			#Extflux2 = Extflux1*Extinc2F
# 
# 		if self._ext_model == 'CCM89':
# 			Extinc1 = self.CCM89(self._wavelengths/1e4, self._extinction_1, self._reddening_1)
# 			Extinc2 = self.CCM89(self._wavelengths/1e4, self._extinction_2, self._reddening_2)
# 
# 			Extinc1F = 2.5**Extinc1 
# 			Extinc2F = 2.5**Extinc2
# 
# 			Extflux1 = finalflux/Extinc1F
# 			Extflux2 = Extflux1/Extinc2F
# 			
# 		elif self._ext_model == 'VCG04':
# 			Extinc1 = self.VCG04(self._wavelengths/1e4, self._extinction_1, self._reddening_1)
# 			Extinc2 = self.VCG04(self._wavelengths/1e4, self._extinction_2, self._reddening_2)
# 
# 			Extinc1F = 2.5**Extinc1 
# 			Extinc2F = 2.5**Extinc2
# 
# 			Extflux1 = finalflux/Extinc1F
# 			Extflux2 = Extflux1/Extinc2F
# 
# 			#Extinc1 = VCG04(self._reddening_1)
# 			#Extinc1F = Extinc1.extinguish(1/self._wavelengths*1e4, Av =  self._extinction_1)
# 			#Extinc2 = VCG04(self._reddening_2)
# 			#Extinc2F = Extinc2.extinguish(1/self._wavelengths*1e4, Av = self._extinction_2)
# 
# 		elif self._ext_model == 'F19':
# 			Extinc1 = F19(self._reddening_1)
# 			Extinc1F = Extinc1.extinguish(1/self._wavelengths*1e4, Av =  self._extinction_1)
# 			Extinc2 = F19(self._reddening_2)
# 			Extinc2F = Extinc2.extinguish(1/self._wavelengths*1e4, Av = self._extinction_2)			
# 
# 			Extflux1 = finalflux*Extinc1F
# 			Extflux2 = Extflux1*Extinc2F
# 
# 		elif self._ext_model == 'G23':
# 			Extinc1 = G23(self._reddening_1)
# 			Extinc1F = Extinc1.extinguish(1/self._wavelengths*1e4, Av =  self._extinction_1)
# 			Extinc2 = G23(self._reddening_2)
# 			Extinc2F = Extinc2.extinguish(1/self._wavelengths*1e4, Av = self._extinction_2)			
# 
# 			Extflux1 = finalflux*Extinc1F
# 			Extflux2 = Extflux1*Extinc2F
# 
# 		Extfactor2 = np.mean(self._flux[norm])/np.mean(Extflux2[norm])
# 
# 		model_spectrum = Extflux2 * Extfactor2
# 
# 		chisq, redchisq = self.chisquared(self._flux, model_spectrum, err = self._flux_err)
# 
# 		if plot == 1:
# 			self.plot_sed(model_spectrum, redchisq)
# 
# 		return chisq, redchisq, model_spectrum
# =============================================================================

	def redshift(self, model_spectrum, plot = 0):
		corr_pix = np.where((self._wavelengths > 2500) & (self._wavelengths < 3160))

		corr_wavelengths, corr_flux, corr_model_spectrum = self._wavelengths[corr_pix], self._flux[corr_pix], model_spectrum[corr_pix]
		
		def gauss(x, a, m, s):
			return a * np.exp( -((x-m)**2) / (2 * s**2) )

		corr = signal.correlate(corr_model_spectrum, corr_flux, mode='full')
		lags = signal.correlation_lags(corr_model_spectrum.size, corr_flux.size, mode='full')
		lag = lags[np.argmax(corr)]

		upper, lower = np.argmax(corr) + 5, np.argmax(corr) - 5
		xup, xlow = lag + 5, lag - 5

		popt, pcov = curve_fit(gauss, lags[lower:upper], corr[lower:upper])
		cerr = np.sqrt(np.diag(pcov))

		factor = 1.5  # Number of Angstroms per pixel

		lamsys = popt[1] * factor
		lamsys_err = cerr[1] * factor

		if plot == 1:
			plt.rc('axes', labelsize=12)
			plt.rc('figure', titlesize=16)
			plt.rc('xtick', labelsize=12)
			plt.rc('ytick', labelsize=12)
			plt.rc('legend', fontsize=12)

			fig, ax = plt.subplots(1, 1, figsize=(6,4))

			x1 = np.linspace(lag - 5, lag + 5, 1000)

			ax.plot(lags[lower:upper]*factor, corr[lower:upper], color = 'k', label='Correlated Wavelength Lags')
			ax.plot(x1*factor, gauss(x1, *popt), color='r', label='Best Fit Gaussian\n$\mu$=' + str(round(lamsys,3)) + '±' + str(round(lamsys_err,3)))
			ax.set_xlabel('Wavelength Lag (Ang.)')
			ax.set_ylabel('Correlation Strength')

			top = corr[lower:upper].max()
			bottom = gauss(x1,*popt).min()			

			ax.set_ylim(top = top+(top-bottom)*0.5, bottom = bottom)
			ax.tick_params(axis='both', direction='in', which='both')
			ax.legend(loc='upper left')
			ax.set_title(self._star + ': Best Fit Wavelength Lags')	
			fig.tight_layout()
			plt.savefig(self._star + "_Redshift.pdf")
			plt.show()

		return lamsys, lamsys_err

	def interpolate(temps, gravs, metals, dimen, out):
		flux_arr, flux_interp = [],[]
		if dimen == 'temp': # dont think this part works
			for temp in temps:
				fluxfile = '/d/tel2/brock/Catalogs/PHOENIX/lte' + temp.zfill(5) + "-" + gravs[0] + metals[0] + ".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
				wavefile="/d/tel2/brock/Catalogs/PHOENIX/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
				flux, header = fits.getdata(fluxfile, header=True)
				wavs, header2 = fits.getdata(wavefile, header=True)
				flux_arr.append(flux)

			flux_arr = np.array(flux_arr)
			temp_vals = np.array(temps, dtype=np.float32) # input metallicities
			#interp = np.linspace(float(metals[0]),float(metals[-1]), 11) # metallicities to  interpolate at
			interp_func = interp1d(temp_vals, flux_arr, axis=0, kind='quadratic')

			flux_interp = interp_func(out)

			for i in range(np.shape(flux_interp)[0]):
				hdu = fits.PrimaryHDU(flux_interp[i,:])
				write_file = "/d/tel2/brock/Catalogs/PHOENIX/lte" + out[i].zfill(5) + "-" + gravs[0] + metals[0] + ".PHOENIX-ACES-AGSS-COND-2011-HiRes_test.fits"
				hdu.writeto(write_file, overwrite=True)

		elif dimen == 'grav':
			for grav in gravs:
				print(grav)
				fluxfile = '/d/tel2/brock/Catalogs/PHOENIX/lte' + temps[0].zfill(5) + "-" + grav + metals[0] + ".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
				wavefile="/d/tel2/brock/Catalogs/PHOENIX/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
				flux, header = fits.getdata(fluxfile, header=True)
				wavs, header2 = fits.getdata(wavefile, header=True)
				flux_arr.append(flux)

			flux_arr = np.array(flux_arr)
			grav_vals = np.array(gravs, dtype=np.float32) # input metallicities
			interp_func = interp1d(grav_vals, flux_arr, axis=0, kind='quadratic')

			flux_interp = interp_func(out)

			for i in range(np.shape(flux_interp)[0]):
				hdu = fits.PrimaryHDU(flux_interp[i,:])
				write_file = "/d/tel2/brock/Catalogs/PHOENIX/lte" + temps[0].zfill(5) + "-" + out[i] + metals[0] + ".PHOENIX-ACES-AGSS-COND-2011-HiRes_test.fits"
				hdu.writeto(write_file, overwrite=True)

		elif dimen == 'metal':
			for metal in metals:
				fluxfile = '/d/tel2/brock/Catalogs/PHOENIX/lte' + temps[0].zfill(5) + "-" + gravs[0] + metal + ".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
				wavefile="/d/tel2/brock/Catalogs/PHOENIX/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
				flux, header = fits.getdata(fluxfile, header=True)
				wavs, header2 = fits.getdata(wavefile, header=True)
				flux_arr.append(flux)

			flux_arr = np.array(flux_arr)
			metal_vals = np.array(metals, dtype=np.float32) # input metallicities
			interp_func = interp1d(metal_vals, flux_arr, axis=0, kind='linear')

			flux_interp = interp_func(out)

			for i in range(np.shape(flux_interp)[0]):
				hdu = fits.PrimaryHDU(flux_interp[i,:])
				write_file = "/d/tel2/brock/Catalogs/PHOENIX/lte" + temps[0].zfill(5) + "-" + gravs[0] + out[i] + ".PHOENIX-ACES-AGSS-COND-2011-HiRes_test.fits"
				hdu.writeto(write_file, overwrite=True)

	def mcmc(self, av1_range, rv1_range, av2_range, rv2_range):
		print(av1_range)
		self.fit_sed()
		
