"""
Script to interface with extinction.reddening class to fit continuum normalized spectrum with normalized PHOENIX spectrum.
"""

import os
os.chdir('/d/users/brock/python/')
import extinction
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.visualization import quantity_support
# BP Neccessary imports.

quantity_support()
# BP Allowing units to plot properly.

red_spec = extinction.reddening()
# BP Creating reddening object.

name = 'pk04_opt_norm'
# BP Defining target name.

wavelength, flux, flux_err = red_spec.read_spectrum(name, '/d/tel2/brock/Data/Extinction/')
wavelength, flux, flux_err = u.Quantity(wavelength, u.AA), u.Quantity(flux, u.erg / (u.s * u.cm**2 * u.AA)), u.Quantity(flux_err, u.erg / (u.s * u.cm**2 * u.AA))
# BP Extracting HST raw flux for given target in the right units.

wavelength += 1*u.AA

ind = (wavelength > 4200*u.AA) & (wavelength < 4900*u.AA)
wavelength, flux, flux_err = wavelength[ind], flux[ind], flux_err[ind]
# BP Chopping off all wavelengths blue of 2640, just before a large drop.

cont_flux = red_spec.continuum_normalize(wavelength, flux, flux_err)
# BP Continuum normalizing the 

cont_factor = cont_flux/flux
cont_flux_err = flux_err * cont_factor
# BP Renormalizing raw HST flux error to normalized spectrum.

temperature = 5750
gravity = 4.134
metallicity = None
smooth_factor = 5
# BP Defining stellar atmosphere model parameters.

# 2750 doesnt seem to change
# 2650 doesnt seem to change
# 2930 constant but lines around it change a lot

model_temperature, model_temp_err, model_gravity, model_grav_err, model_metallicity, model_metal_err, final_model_spectrum, model_spectrum_lower, model_spectrum_upper = red_spec.fit_stellar_lsq(wavelength, cont_flux, cont_flux_err, smooth_factor, temp_in = temperature, grav_in = gravity, metal_in = metallicity)
# BP Fitting stellar spectrum with PHOENIX model using least squares method.
# BP Any stellar parameters (temperature, metallicidty, gravity) not specified or given as None will be fit.
# BP Pass specific stellar parameters to hold them fixed.

chisq, redchi = red_spec.chi_squared(cont_flux, final_model_spectrum, cont_flux_err)
# BP Calculating chi squared of the model.

fig, ax = plt.subplots(figsize=(10,4), layout='tight')
# BP Creating plotting axes.

ax.plot(wavelength, cont_flux, label = 'HST', color='b')
ax.plot(wavelength, final_model_spectrum, label = '$\chi^2_{{red}}$: {:.3f}'.format(redchi), color= 'r')
ax.fill_between(wavelength, model_spectrum_lower, model_spectrum_upper, facecolor = 'r', alpha = 0.35)
# BP Plotting continuum normalized HST spectrum and best fit model.

ax.axhline(1, color = 'k')

ax.legend(loc = 'lower left')
ax.text(4900, 0.1, '{}:\n$T_{{eff}}$={:.1f}, [Fe/H]=${:.3f}$, log $g$={:.3f}\n±{:.1f},             ±{:.3f},        ±{:.3f}'.format(name, model_temperature, model_metallicity, model_gravity, model_temp_err, model_metal_err, model_grav_err), size = 18, horizontalalignment='right')
ax.set_ylabel('Normalized Flux')
ax.set_ylim([0, 1.2])
plt.savefig('/d/users/brock/paper_figures/{}_stellar.pdf'.format(name))
plt.show()
# BP Displaying best fit parameters and axes.
