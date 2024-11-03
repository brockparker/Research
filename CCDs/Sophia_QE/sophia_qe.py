#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sophia QE analysis script. Combines functionality from MC_characterization_script.py, qe_analysis_v2.py, and qe_photodiode_analysis.py from Aafaque Khan. 
Simplified data reduction for testing purposes.

@author: baparker
"""

import re
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.units as u
import astropy.constants as c
from astropy.visualization import quantity_support
import scipy.interpolate
from astropy.io import fits
from matplotlib import colormaps
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

path_base = '/home/baparker/GitHub/'

data_path = path_base + 'Research/CCDs/Sophia_QE/20240607/'
# BP Directory storing raw data, both fits files from sophia (dark, bias, and science) and photodiode csvs (dark and science).
photodiode_response_file = path_base + 'Research/CCDs/Sophia_QE/photodiode_response_new.csv'
# BP File location of the csv file containing the current to power conversion for the calibrated photodiode.
# BP This will vary for different photodiodes.
compiled_file_path = data_path + 'compiled_data.csv'
# BP File location of the master database of all compiled photodiode and fits file region data. All plots and analysis can be done from this file.

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
    # BP Get the unsigned photodiode current.
    
    response = get_photodiode_response(photodiode_response_file)
    conversion = response(wavelength.to(u.m).value) * u.A / u.W
    # BP Calculate the photodiode conversion factor at each wavelength.
    
    photodiode_power = current / conversion
    # BP Convert photodiode current into power.

    photon_energy = c.h * c.c / wavelength    
    # BP Calculate the energy of a photon at each wavelength.

    photon_rate = photodiode_power / photon_energy * u.photon
    # BP Calculate the rate of photons needed to produce the observed power.

    return photon_rate.to(u.photon / u.s)

def compile_photodiode_csv(data_path, compiled_data_path, string_pattern, dark_string_pattern = 0):
    
    df = pd.DataFrame(data=None, columns=['wavelength', 'ch1_mean', 'photon_rate', 'sophia_region'], dtype='float64')
    
    photodiode_in_list = glob.glob(data_path + string_pattern)
    
    if dark_string_pattern != 0:
        dark_in_list = glob.glob(data_path + dark_string_pattern)
        dark_data = pd.read_csv(dark_in_list)
        
        photodiode_dark = np.mean(dark_data['Ch1'])
        
    else:
        photodiode_dark = 0
    
    for idx, file in enumerate(photodiode_in_list):
        data = pd.read_csv(file)
        
        ind = re.search(string_pattern.replace('.','\.').replace('???','...'), file)
        start = ind.start() + 6
        end = ind.end() - 4
        
        wavelength = int(file[start:end])
        # BP Wavelength in nanometers
        
        ch1 = data['Ch1']
        ch1_mean = np.mean(ch1) - photodiode_dark
        
        compiled_data = np.array([wavelength, ch1_mean])
        compiled_data = np.sort(compiled_data)
        
        df.loc[idx, ['wavelength', 'ch1_mean']] = np.array([wavelength, ch1_mean])

    final_data = update_master_database(compiled_file_path, df)
    
    return final_data

def reduce_camera_images_average(data_path, bias_string_pattern, dark_string_pattern, science_string_pattern):
    bias_in_list = glob.glob(data_path + bias_string_pattern)
    dark_in_list = glob.glob(data_path + dark_string_pattern)
    science_in_list = glob.glob(data_path + science_string_pattern)
    bias_in_list.sort()
    dark_in_list.sort()
    science_in_list.sort()
    
    # Should probably read in and store these all as a 3D datacube.
    master_bias, master_dark = 0, 0
    
    # BP Create master bias.
    for file in bias_in_list:
        hdul = fits.open(file)
        hdr = hdul[0].header
        data = hdul[0].data[0]
        
        master_bias += data
        
        hdul.close()
        
    master_bias = master_bias / len(bias_in_list)
    
    hdu = fits.PrimaryHDU(data = master_bias, header = hdr)
    hdu.writeto(data_path + 'master_bias.fits', overwrite = True)
    
    dark_bias_subtracted_list = []
    
    # BP Subtract master bias from all dark images.
    for file in dark_in_list:
        outfile = file.replace('.fits', '_b.fits')
        
        hdul = fits.open(file)
        hdr = hdul[0].header
        data = hdul[0].data[0]
        
        data = data - master_bias
        
        hdu = fits.PrimaryHDU(data = data, header = hdr)
        hdu.writeto(outfile, overwrite = True)
        
        dark_bias_subtracted_list.append(outfile)
        
        hdul.close()
        
    science_bias_subtracted_list = []

    # Subtract master bias from all science images.
    for file in science_in_list:
        outfile = file.replace('.fits', '_b.fits')
        
        hdul = fits.open(file)
        hdr = hdul[0].header
        data = hdul[0].data[0]
        
        data = data - master_bias
        
        hdu = fits.PrimaryHDU(data = data, header = hdr)
        hdu.writeto(outfile, overwrite = True)
        
        science_bias_subtracted_list.append(outfile)
        
        hdul.close()
    
    # BP Create master dark by scaling all darks to 1 second after subtracting bias.
    for file in dark_bias_subtracted_list:       
        hdul = fits.open(file)
        hdr = hdul[0].header
        data = hdul[0].data[0]
        
        exp_time = float(hdr['EXPTIME'])
        
        try:
            scaled_data = data / exp_time
        except:
            print('Temporary issue with dark measurements')
            scaled_data = data
        
        master_dark += scaled_data
        
        hdul.close()
        
    master_dark = master_dark / len(dark_bias_subtracted_list)
        
    hdu = fits.PrimaryHDU(data = master_dark, header = hdr)
    hdu.writeto(data_path + 'master_dark.fits', overwrite = True)
    # BP Currently generating and using master_dark, may be better to use individual darks for each wavelength, as exposure time should be the same.
    
    science_out_list = []
    
    # BP Subtract master dark from all science images
    for file in science_bias_subtracted_list:
        outfile = file.replace('_b.fits', '_bd.fits')
        
        hdul = fits.open(file)
        hdr = hdul[0].header
        data = hdul[0].data[0]
        
        data = data - master_dark
        
        hdu = fits.PrimaryHDU(data = data, header = hdr)
        hdu.writeto(outfile, overwrite = True)
        
        science_out_list.append(outfile)
        
        hdul.close()
    
    return science_out_list

def reduce_camera_images_individual(data_path, bias_string_pattern, dark_string_pattern, science_string_pattern):
    bias_in_list = glob.glob(data_path + bias_string_pattern)
    dark_in_list = glob.glob(data_path + dark_string_pattern)
    science_in_list = glob.glob(data_path + science_string_pattern)
    
    bias_in_list.sort()
    dark_in_list.sort()
    science_in_list.sort()
    
    for b_file, d_file, s_file in zip(bias_in_list, dark_in_list, science_in_list):
        hdul_b = fits.open(b_file)
        data_b = hdul_b[0].data
        
        hdul_d = fits.open(d_file)
        data_d = hdul_d[0].data
        
        hdul_s = fits.open(s_file)
        hdr_s = hdul_s[0].header
        data_s = hdul_s[0].data
        
        science_outfile = s_file.replace('.fits', '_bd_separate.fits')
        
        dark_sub_bias = data_d - data_b
        science_sub_bias = data_s - data_b
        
        science_out = science_sub_bias - dark_sub_bias
        
        hdu = fits.PrimaryHDU(data = science_out, header = hdr_s)
        hdu.writeto(science_outfile, overwrite = True)
    
def update_master_database(compiled_file_path, df):
    # Helper function to update a specific row in the master compiled data csv.
    print('Updating compiled data.')
# =============================================================================
#     if os.path.isfile(compiled_file_path):
#         final_data = pd.read_csv(compiled_file_path)
#         
#         return final_data
#             
#     else:
# =============================================================================
    df = df.sort_values('wavelength', ignore_index = True)
    df.to_csv(compiled_file_path, index = False)
    
    return df

def main():   
    print('Reducing FITS files.')
    sophia_file_list = reduce_camera_images_individual(data_path, 'bias_???_2.fits', 'dark_???_2.fits', 'science_???_2.fits')
    ### Reduce raw sophia files. Subtract biases, subtract (and scale if necessary) darks.
    
    print('Compiling photodiode data.')
    ### Compile all raw photodiode csvs into one master photodiode csv. Take the mean of each raw file.
    global photodiode_photon_rate, wavelengths, compiled_data_frame, photon_rate_df
    compiled_data_frame = compile_photodiode_csv(data_path, compiled_file_path, 'picoa_???.csv')
    print(compiled_data_frame)


    print('Converting photodiode response into photon counts.')
    # BP Read in photodiode response function.
    ### Convert master photodiode csv of currents into total photon counts. Read in conversion and do calculations.
    wavelengths = np.array(compiled_data_frame['wavelength']) * u.nm
    photodiode_current = np.array(compiled_data_frame['ch1_mean']) * u.A
    
    photodiode_photon_rate = current_to_photon_rate(wavelengths, photodiode_current)
    
    photon_rate_df = pd.DataFrame(data = photodiode_photon_rate, columns = ['new'], dtype='float64')


    compiled_data_frame.update(photon_rate_df)
    print(compiled_data_frame)
    #compiled_data_frame.loc[wavelengths.value == compiled_data_frame['wavelength'], ['photon_rate']] = np.array([photodiode_photon_rate])
    
    print('Reading in sophia regions.')
    ### Add mean counts of different regions in sophia data into new file with wavelengths and photodiode values.
    
if __name__ == "__main__":
    main()

