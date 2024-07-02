#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sophia QE analysis script. Combines functionality from MC_characterization_script.py, qe_analysis_v2.py, and qe_photodiode_analysis.py from Aafaque Khan. 
Simplified data reduction for testing purposes.

@author: baparker
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.units as u
import astropy.constants as c
from astropy.visualization import quantity_support
import scipy.interpolate
# BP Necessary imports.

quantity_support()

plt.rc('axes', labelsize=14)
plt.rc('figure', titlesize=30)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('legend', fontsize=12)
plt.rc('xtick.major', size=6)    # size of the tick markers
plt.rc('ytick.major', size=6)    # size of the tick markers
plt.rc('xtick.minor', size=4)    # size of the tick markers
plt.rc('ytick.minor', size=4)    # size of the tick markers
# BP Plot stylization parameters.

data_path = '/home/baparker/GitHub/Research/CCDs/Sophia_QE/20240607/'
# BP Directory storing raw data, both fits files from sophia (dark, bias, and science) and photodiode csvs (dark and science).
photodiode_response_file = '/home/baparker/GitHub/Research/CCDs/Sophia_QE/photodiode_response_new.csv'
# BP File location of the csv file containing the current to power conversion for the calibrated photodiode.
# BP This will vary for different photodiodes.

def log_interp1d(xx, yy, kind='linear'):
    """
    Helper function to do logarithmic interpolation. Converts into logarithm space, does the interpolation, then converts back into linear space.

    Parameters
    ----------
    xx : np.array
        Input array of x values to interpolate.
    yy : np.array
        Input array of y values to interpolate.
    kind : string, optional
        The type of interpolation to do. The default is 'linear'.

    Returns
    -------
    log_spl : scicpy.interpolate._fitpack2.LSQUnivariateSpline
        Interpolated photodiode response function. Input wavelength to get out the corresponding conversion factor.

    """
    logx = np.log10(xx)
    logy = np.log10(yy)
    # BP Take the logarithm of both variables.
    
    lin_spl = scipy.interpolate.UnivariateSpline(logx, logy)
    lin_spl.set_smoothing_factor(0.0075)
    log_spl = lambda zz: np.power(10.0, lin_spl(np.log10(zz)))
    # BP Take the interpolation in logarithmic space, then convert back into linear space.
    
    # lin_interp = scipy.interpolate.interp1d(logx, logy, kind=kind)
    # log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
    # BP Unused interpolation function.
    return log_spl

def get_photodiode_response(photodiode_response_file):
    """
    Function to get the photodiode current to power conversion at a given wavelength and interpolate.
    The interpolation is only done once and recalled if available.

    Parameters
    ----------
    wavelength : astropy.Quantity
        Array-like astropy quantity of input wavelength. Must be in units convertible to meters.
    photodiode_response : string
        File location of the file containing the photodiode response function.

    Returns
    -------
    response_interpolated : scicpy.interpolate._fitpack2.LSQUnivariateSpline
        Interpolated photodiode response function. Input wavelength to get out the corresponding conversion factor.

    """
    
    response = pd.read_csv(photodiode_response_file)
    # BP Read in CSV containing current to power conversion for the photodiode.
    
    x = np.array(response['wavelength'])
    # BP Get wavelength array, assuming it is in meters.
    
    y = np.array(response['conversion'])
    # BP Get photodiode conversion factor. Should be a value between 0 and 1.
    
    response_interpolated = log_interp1d(x, y, kind='quadratic')
    # BP Interpolate the loaded x and y values. Done in logarithmic space and then converted back into linear space.
    
    return response_interpolated

@u.quantity_input(wavelength=u.m, current=u.A)
def current_to_photon_rate(wavelength, current):
    """
    Function to convert the measured current from the photodiode to photon count rate.

    Parameters
    ----------
    wavelength : astropy.Quantity
        Wavelength of measured photocurrent. Must be in units convertible to meters.
    current : astropy.Quantity
        Signed photocurrent measured from the VUV photodiode. Must be in units convertible to amps.

    Returns
    -------
    photon_rate : astropy.Quantity
        Returns the photon rate as an astropy Quantity.
    """
    
    current = np.abs(current)
    # BP Get the unsigned current.
    
    response = get_photodiode_response(photodiode_response_file)
    conversion = response(wavelength.to(u.m).value) * u.A / u.W
    
    photodiode_power = current / conversion

    photon_energy = c.h * c.c / wavelength    

    photon_rate = photodiode_power / photon_energy * u.photon

    return photon_rate.to(u.photon / u.s)

def main():


    
    print('Reading in raw data.')
    ### Read in raw data
    
    print('Reducing FITS files.')
    ### Reduce raw sophia files. Subtract biases, subtract (and scale if necessary) darks.
    
    print('Compiling photodiode data.')
    ### Compile all raw photodiode csvs into one master photodiode csv. Take the mean of each raw file.
    
    #current = current * u.A
    
    print('Converting photodiode response into photon counts.')
    # BP Read in photodiode response function.
    ### Convert master photodiode csv of currents into total photon counts. Read in conversion and do calculations.
    
    print('Reading in sophia regions.')
    ### Add mean counts of different regions in sophia data into new file with wavelengths and photodiode values.
    
if __name__ == "__main__":
    main()

