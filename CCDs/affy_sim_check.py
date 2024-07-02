#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 10:48:35 2024

@author: baparker
"""

import numpy as np
import astropy.constants as c
import astropy.units as u
from astropy.modeling import models
import matplotlib.pyplot as plt
from astropy.visualization import quantity_support

quantity_support()

def get_bb_flux_per_pix(T=298*u.K, pix_size=13e-6*u.m, wlb=200e-9*u.m, wlr=1000e-9*u.m):
    wavelengths = np.linspace(wlb, wlr, 1000)  # Wavelengths from 200 nm to 800 nm
    area = pix_size**2  # Area in square meters


    planck_me = models.BlackBody(temperature = T)

    def planck(wavelength, T):
        return (2 * c.h * c.c**2) / (wavelength**5) / (np.exp((c.h * c.c) / (wavelength * c.k_B * T)) - 1)

    spectral_radiance = planck(wavelengths, T)
    spectral_radiance_photons = spectral_radiance * wavelengths / (c.h * c.c)    
    flux_photons_per_sr = np.trapz(spectral_radiance_photons, wavelengths)
    flux_photons = flux_photons_per_sr * 2 * np.pi
    total_flux_photons = flux_photons * area
    
    return total_flux_photons


wavs = np.linspace(200, 100000, 1000)*u.nm
T = 298 * u.K
planck_me = models.BlackBody(temperature = T)

def planck(wavelength, T):
    #return (2 * c.h * c.c**2) / (wavelength**5) / (np.exp((c.h * c.c) / (wavelength * c.k_B * T)) - 1)
    return 2 * c.h * c.c**2 / (wavelength**5 * (np.exp(c.h * c.c / (wavelength * c.k_B * T)) - 1))

plot_me = planck_me(wavs)
plot_affy = planck(wavs, T)

plt.plot(wavs, plot_me.si, label='Me')
plt.plot(wavs, plot_affy.si, label='Affy')
plt.legend(loc = 'best')
plt.show()