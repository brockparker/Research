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
# BP Necessary imports.

quantity_support()
# BP Enable astropy quantity support.

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
###path_base = r'C:\Users\Brock\Documents\Git'

data_path = path_base + 'Research/CCDs/Sophia/Data/20241102/QE/'
###data_path = path_base + r'\Research\CCDs\Sophia\Data\20240607\\'
# BP Directory storing raw data, both fits files from sophia (dark, bias, and science) and photodiode csvs (dark and science).
photodiode_response_file = path_base + 'Research/CCDs/Sophia/QE/photodiode_response_new.csv'
###photodiode_response_file = r'C:\Users\Brock\Documents\Git\Research\CCDs\Sophia\QE\photodiode_response_new.csv'
# BP File location of the csv file containing the current to power conversion for the calibrated photodiode.
# BP This will vary for different photodiodes.
#port_relation_file = N/A
# BP File location of the relation between main port and side port on the monochromator.
# BP This should be standard for the monochromator setup, although may change will distance of detector from the port.
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
    photodiode_response_file : string
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

def get_port_relation(port_relation_file = None):
    """
    Function to get the correction factor between the main port and side port of the monochromator, used to correct the photodiode counts at the side port.
    The interpolation is only done once and recalled if available.

    Parameters
    ----------
    port_relation_file : string
        File location of the file containing the port relation file.

    Returns
    -------
    response_interpolated : scicpy.interpolate._fitpack2.LSQUnivariateSpline
        Interpolated photodiode response function. Input wavelength to get out the corresponding port relation.
    """
    if port_relation_file:        
        ##############################################################
        # TODO: Will need to update when the actual data is recovered.
        ##############################################################
        
        relation = pd.read_csv(port_relation_file)
        # BP Read in CSV containing side port to main port relation.
        
        x = np.array(relation['wavelength'])
        # BP Get wavelength array, assuming it is in meters.
        
        y = np.array(relation['relation'])
        # BP Get photodiode conversion factor. Should be a value between 0 and 1.
        
        relation_interpolated = log_interp1d(x, y, kind='quadratic')
        # BP Interpolate the loaded x and y values. Done in logarithmic space and then converted back into linear space.
        
        return relation_interpolated
    
    else:
        relation_interpolated = lambda x: 1
        # BP Simply return a factor of unity, i.e. no correction, if no file is specified.
        
        return relation_interpolated

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
        
    relation = get_port_relation()
    port_conversion = relation(wavelength.to(u.m).value)
    # BP Calculate port relation correction factor at each wavelength.
        
    corrected_current = current / port_conversion
    # BP Convert into main port flux from side port flux.
    
    photodiode_power = corrected_current / conversion
    # BP Convert photodiode current into power.

    photon_energy = c.h * c.c / wavelength    
    # BP Calculate the energy of a photon at each wavelength.

    photon_rate = photodiode_power / photon_energy * u.photon
    # BP Calculate the rate of photons needed to produce the observed power.

    return photon_rate.to(u.photon / u.s)

def compile_photodiode_csv(data_path, compiled_data_path, string_pattern, dark_string_pattern = None):
    """
    Parameters
    ----------
    data_path : string
        File path to the folder containing all photodiode data.
    compiled_data_path : string
        File path to location where compiled data frame csv is stored. Typically inside of data path.
    string_pattern : string
        Unix wildcard style string to search for photodiode csv files.
        Standard is 'picoa_???.csv'
    dark_string_pattern : string, optional
        Unix wildcard style string to search for photodiode csv files.
        Standard is 'picoa_0.csv'
        
    Returns
    -------
    final_data : pd.dataFrame
        Pandas data frame containting the final compiled data from the photodiode and camera regions. This is updated every step of the reduction.
    """
    
    df = pd.DataFrame(data=None, columns=['wavelength', 'ch1_mean', 'ch1_std', 'photon_rate', 'photon_rate_std', 'sophia_region', 'sophia_region_std', 'qe'], dtype='float64')
    # BP Create blank dataframe to store data.
    
    photodiode_in_list = glob.glob(data_path + string_pattern)
    # BP Find all picoameter csv files.
    
    if dark_string_pattern:
        dark_in_list = glob.glob(data_path + dark_string_pattern)
        dark_data = pd.read_csv(dark_in_list)
        
        photodiode_dark = np.mean(dark_data['Ch1'])
        photodiode_dark_std = np.std(dark_data['Ch1'])
        # BP If a photodiode dark is provided, record its mean and std.
        
    else:
        photodiode_dark, photodiode_dark_std = 0, 0
        # BP Otherwise, set the dark rate to 0 since it should be negligible.
    
    for idx, file in enumerate(photodiode_in_list):
        data = pd.read_csv(file)
        # BP Read in data from photodiode csv.
        
        ind = re.search(string_pattern.replace('.','\.').replace('???','...'), file)
        start = ind.start() + 6
        end = ind.end() - 4
        # BP Look for a specific pattern in the 
        
        wavelength = int(file[start:end])
        # BP Extract wavelength in nanometers from photodiode file name.
                
        ch1_mean = np.mean(data['Ch1'])
        ch1_true = ch1_mean - photodiode_dark
        ch1_std = np.std(data['Ch1']) - photodiode_dark_std
        # TODO: Is this correct
        # BP Photodiode data is the mean of the first channel, which is the only one plugged in. Subtract from it the dark/bias level which should be a separate file.
        
        compiled_data = np.array([wavelength, ch1_true, ch1_std])
        # BP Compile data into array.
                
        df.loc[idx, ['wavelength', 'ch1_mean', 'ch1_std']] = compiled_data
        # BP Update dataframe with wavelength and picoammeter data.
                    
    final_data = update_master_database(df, list(df.columns), compiled_file_path)
    # BP Update master dataframe rows of all input wavelengths.
        
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
        data = hdul[0].data
        
        exp_time = float(hdr['EXPTIME'])
        
        with np.errstate(divide='raise'):  
            try:
                scaled_data = data / exp_time
                master_dark += scaled_data

            except FloatingPointError:
                print('Divide by zero error encountered in exposure time. Double check dark files!')
                master_dark = np.zeros(np.shape(data))
            except:
                print('Temporary issue with dark measurements')
                master_dark = np.zeros(np.shape(data))
            
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
        data = hdul[0].data
            
        exp_time = float(hdr['EXPTIME'])
        
        data = data - master_dark * exp_time
        # BP Scale master dark by exposure time of the science image.
        
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
        
        ### TODO: Need to save these files to read in 
    
def extract_regions(wavelengths, file_list, regions, compiled_file_path):
    df = pd.DataFrame(data=None, columns=['wavelength', 'sophia_region', 'sophia_region_std'], dtype='float64')
    # BP Create blank dataframe to store data.
    
    for idx, file in enumerate(file_list):
        hdul = fits.open(file)
        data = hdul[0].data
        
        wavelength = wavelengths[idx]
        
        for region in regions:            
            region_data = data[region[0]:region[1], region[2]:region[3]]
            
            region_mean = np.mean(region_data)
            region_std = np.std(region_data)
        
            compiled_data = np.array([wavelength.value, region_mean, region_std])

            df.loc[idx, ['wavelength', 'sophia_region', 'sophia_region_std']] = compiled_data
            # BP Update dataframe with wavelength and picoammeter data.
                        
    final_data = update_master_database(df, list(df.columns), compiled_file_path)
    # BP Update master dataframe rows of all input wavelengths.
            
    return final_data

def calculate_qe(compiled_df, compiled_file_path):
    photon_rate = compiled_df['photon_rate']
    photon_rate_stde = compiled_df['photon_rate_std']
    sophia_rate = compiled_df['sophia_region']
    sophia_rate_std = compiled_df['sophia_region_std']
    
    qe = sophia_rate/photon_rate
    
    final_data = update_master_database(qe, ['qe'], compiled_file_path)
    # BP Update master dataframe rows of all input wavelengths.
    
    return final_data
    
def update_master_database(data, columns, compiled_file_path):
    """
    Parameters
    ----------
    data : np.array OR pd.dataFrame
        Array-like set of input data to append to the existing compiled data frame.
    columns : list
        List of strings containing columns to append new data to, e.g. ['wavelength', 'photon_rate']
    compiled_file_path : string
        DESCRIPTION.

    Returns
    -------
    master_df : TYPE
        DESCRIPTION.
    """
    
    # Helper function to update a specific row in the master compiled data csv.
    print('Updating compiled data.')
        
    new_df = pd.DataFrame(data = data, columns = columns, dtype='float64')
    # BP Calculate corresponding photon flux and save it into new data frame to add onto total compiled data frame.

    if os.path.exists(compiled_file_path):
        master_df = pd.read_csv(compiled_file_path)

        master_df.update(new_df)
        # BP Update the compiled data frame with new photon rates.
    else:
        master_df = new_df
        # BP Save blank data frame csv.

    master_df = master_df.sort_values('wavelength', ignore_index = True)
    master_df.to_csv(compiled_file_path, index = False)
    ### TODO TOtally borken
    return master_df

def main():   
# =============================================================================
#     print('Reducing FITS files.')
#     global compiled_data_frame, sophia_file_list
# 
#     sophia_file_list = reduce_camera_images_average(data_path, 'bias_???.fits', 'dark_???.fits', 'science_???_2.fits')
#     ### Reduce raw sophia files. Subtract biases, subtract (and scale if necessary) darks and return list of reduced science files.
#     
#     print('Compiling photodiode data.')
#     compiled_data_frame = compile_photodiode_csv(data_path, compiled_file_path, 'picoa_???.csv')
#     ### Compile all raw photodiode csvs into one master photodiode csv. Take the mean and std of each raw file.
# 
#     print('Converting photodiode response into photon counts.')
#     ### Convert master photodiode csv of currents into total photon counts. Read in conversion and do calculations.
#     wavelengths = np.array(compiled_data_frame['wavelength']) * u.nm
#     photodiode_current = np.array(compiled_data_frame['ch1_mean']) * u.A
#     photodiode_current_std = np.array(compiled_data_frame['ch1_std']) * u.A
#     # BP Extract wavelengths and currents from compiled photodiode csv.
#         
#     photon_rate = current_to_photon_rate(wavelengths, photodiode_current)    
#     photon_rate_std = current_to_photon_rate(wavelengths, photodiode_current_std)
#     # BP Calculate photon_rate from photodiode currents.
#     
#     compiled_data_frame = update_master_database(np.transpose([wavelengths, photon_rate, photon_rate_std]), ['wavelength', 'photon_rate', 'photon_rate_std'], compiled_file_path)
#     # BP Update the master compiled data csv with photon rates.
# =============================================================================
    
    compiled_file_path = '/home/baparker/GitHub/Research/CCDs/Sophia/Data/QE/compiled_data.csv'

    compiled_data_frame = pd.read_csv(compiled_file_path)

    global sophia_file_list

    sophia_file_list = glob.glob(path_base + 'Research/CCDs/Sophia/Data/QE/*.fits')
    sophia_file_list.sort()
    wavelengths = np.array(compiled_data_frame['wavelength']) * u.nm
    
    print(sophia_file_list)
    
    print('Reading in sophia regions.')
    region_array = np.array([[[0, 2048, 0, 2048]]
                        , [[500, 750, 500, 750]]
                        , [[1000, 1250, 1000, 1250]]
                        , [[1250, 1750, 1250, 1750]]])
                        
    
    for regions in region_array:
        print(regions)
        compiled_data_frame = extract_regions(wavelengths, sophia_file_list, regions, compiled_file_path)
        ### Add mean counts of different regions in sophia data into new file with wavelengths and photodiode values.
        
        print('Calculating Quantum Efficiency')
        compiled_data_frame = calculate_qe(compiled_data_frame, compiled_file_path)
        
        qe = compiled_data_frame['qe']
        plt.plot(wavelengths, qe)
        plt.xlim(250*u.nm, 1050*u.nm)
        plt.show()
        
if __name__ == "__main__":
    main()

