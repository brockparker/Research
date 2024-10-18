import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter

######################################################################################################
lam1, lam2, cdelt1 = 	1400, 25000, 0.1		# wavelength range
sigma = 		20				# amount of smoothing
red = 			3.1
temp = 			'10000'
grav = 			'4.00'
ext = 			1.00
metal = 		'-0.0'
######################################################################################################
def CCM89(microns, Av, Rv):
	inverse = np.array(1./microns)
	IR = np.array(np.nonzero(inverse < 1.1)).ravel()
	Optical = np.array(np.nonzero((inverse >= 1.1) & (inverse <= 3.3))).ravel()
	UV = np.array(np.nonzero((inverse > 3.3) & (inverse <= 8.))).ravel()

	if len(IR) > 0:
		a_IR = 0.574 * (inverse[IR]**1.61)
		b_IR = -0.527 * (inverse[IR]**1.61)
		lam_IR = (a_IR + b_IR/Rv) * Av
	else:
		lam_IR = np.array([])

	if len(Optical) > 0:
		y = inverse[Optical]-1.82
		a_opt = 1. + 0.17699*y - 0.50447*(y**2) - 0.02427*(y**3) + 0.72085*(y**4) + 0.01979*(y**5) - 0.77530*(y**6) + 0.32999*(y**7)
		b_opt = 1.41338*y + 2.28305*(y**2) + 1.07233*(y**3) - 5.38434*(y**4) - 0.62251*(y**5) + 5.30260*(y**6) - 2.09002*(y**7)
		lam_opt=(a_opt + b_opt/Rv) * Av
	else:
		lam_opt=np.array([])
	
	if len(UV) > 0:
		UV_2 = np.array(np.nonzero((inverse >= 5.9) & (inverse <= 8.))).ravel()
		UV_1 = np.array(np.nonzero((inverse < 5.9) & (inverse > 3.3))).ravel()

		Fa_2 = -0.04473*(inverse[UV_2]-5.9)**2 - 0.009779*(inverse[UV_2]-5.9)**3
		Fb_2 = 0.2130*(inverse[UV_2]-5.9)**2 + 0.1207*(inverse[UV_2]-5.9)**3
		
		Fa_1 = UV_1 - UV_1	
		Fb_1 = Fa_1		

		F_a = np.concatenate([Fa_2, Fa_1])
		F_b = np.concatenate([Fb_2, Fb_1])

		a_UV = 1.752 - 0.316*inverse[UV] - 0.104/((inverse[UV]-4.67)**2 + 0.341) + F_a
		b_UV = -3.090 + 1.825*inverse[UV] + 1.206/((inverse[UV]-4.62)**2 + 0.263) + F_b
   
		lam_UV = (a_UV + b_UV/Rv) * Av
	else:
		lam_UV=np.array([])

	A_lam = np.concatenate([lam_UV,lam_opt,lam_IR])

	return A_lam

plt.rc('axes', labelsize=16)
plt.rc('figure', titlesize=16)
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('legend', fontsize=16)

fig, ax = plt.subplots(figsize=(8,6), layout='tight')

wavefile="/d/tel2/brock/Catalogs/PHOENIX/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"

fluxfile = '/d/tel2/brock/Catalogs/PHOENIX/lte' + temp.zfill(5) + "-" + grav + metal + ".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
flux, header = fits.getdata(fluxfile, header=True) 	# fluxes in erg/s/cm^2/cm
wavs, header2 = fits.getdata(wavefile, header=True) 	# wavelengths in Ang

flux = flux*1e-8

finalwavs = np.arange(lam1,lam2,cdelt1)	

cs = CubicSpline(wavs, flux)
newflux = cs(finalwavs)
finalflux = gaussian_filter(newflux, sigma)

Extinc1 = CCM89(finalwavs/1e4, ext, red) 	# extinct in mags
Extinc1F = 2.5**Extinc1 				# factor in linear units

Extflux = finalflux/Extinc1F 
spectrum = Extflux

#plt.plot(np.log10(finalwavs), np.log10(finalflux), color = 'k', label = 'Raw Spectrum')
#plt.plot(np.log10(finalwavs), np.log10(spectrum), color = 'r', label = 'Reddened Spectrum')

ax.plot(np.log10(finalwavs), finalflux, color = 'k', label = 'Raw Spectrum')
ax.plot(np.log10(finalwavs), spectrum, color = 'r', label = 'Reddened Spectrum')

ax.text(3.80, 0.61e8, '$T_{{eff}}$=10000, [Fe/H]=0.0,\nlog $g$=4.00', size = 18)
ax.text(3.80, 0.50e8, 'A$_{{V}}$=1.00, R$_{{V}}$=3.1', size = 18, color='r')

ax2 = ax.twiny()
tick_loc = [4.0, 3.65, 3.25]
tick_labels = ['IR', 'Visible', 'UV']
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(tick_loc)
ax2.set_xticklabels(tick_labels)

ax.set_xlabel('log Wavelength (Ang.)')
ax.tick_params(axis='both', direction='in', which='both')
ax.legend(loc='best')
ax.set_ylabel('F$_{\\lambda}$ (erg s$^{-1}$ cm$^{-2}$ Ang$^{-1}$) ')
#ax.set_title('Extincted Stellar Spectrum')
plt.tight_layout()
plt.show()


