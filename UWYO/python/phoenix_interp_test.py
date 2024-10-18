import os
os.chdir('/d/users/brock/python/')
import extinction
import numpy as np
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

plt.rc('axes', labelsize=16)
plt.rc('figure', titlesize=16)
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('legend', fontsize=16)

fig, ax = plt.subplots(figsize=(6,4), layout='tight')

red_spec = extinction.reddening()

wavefile = "/d/tel2/brock/Catalogs/PHOENIX/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
wavelength, header2 = fits.getdata(wavefile, header=True)

ind = (wavelength > 4250) & (wavelength < 4400)
wavelength = wavelength[ind]

temperature = 6000
gravity = 4.00
smooth_factor = 5

for metallicity in np.arange(-0.5, 2.1, 0.5):
	metallicity *= -1
	cs_smooth = red_spec.get_phoenix_sed(temperature, gravity, metallicity, smooth_factor)
	
	flux = cs_smooth(wavelength)
		
	#cont_flux = red_spec.continuum_normalize(wavelength*u.AA, flux * u.erg / (u.s * u.cm**2 * u.AA), np.repeat(flux.mean()*0.15, len(flux)))
		
	ax.plot(wavelength, np.log10(flux), label = metallicity)

ax.set_title('Metallicity Interpolation')
ax.legend(loc = 'best')
plt.margins(x = 0)
plt.show()

fig, ax = plt.subplots(figsize=(6,4), layout='tight')

gravity = 4.00
metallicity = -0.0
smooth_factor = 20

for temperature in np.arange(5500, 7050, 250):
	cs_smooth = red_spec.get_phoenix_sed(temperature, gravity, metallicity, smooth_factor)
	
	flux = cs_smooth(wavelength)
				
	ax.plot(wavelength, np.log10(flux), label = temperature)

ax.set_title('Temperature Interpolation')
ax.legend(loc = 'best')
plt.margins(x = 0)
plt.show()