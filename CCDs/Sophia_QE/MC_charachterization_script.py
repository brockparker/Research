"""
MC_charachterisation_analysis - Python script for analysis of photocurrent measurements for photodidoes mounted on the two ports of the UV-VIS Detector Characterisation Setup at Hamden UV/Vis Detector lab. This scripts generates a database from the simultaneous measurements at the two ports. It converests the photocurrents to photon counts based on the wavelenght of the scan and stores the photocurrent and statistics of the measurements. The ratio between the photocurrent flux at the two ports of the  setup are stored for data reduction for QE measurements. 
"""

__version__ = "0.0.2"

import numpy as np

#import sys
#import warnings
import os
import glob
#import datetime
#from IPython.display import HTML

import matplotlib.pyplot as plt
#from matplotlib.pyplot import cm
#from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

import pandas as pd

#from astropy.visualization import SqrtStretch
#from astropy.visualization.mpl_normalize import ImageNormalize
#from matplotlib.colors import LogNorm
#from astropy.stats import sigma_clipped_stats
#from astropy.modeling import models, fitting
from scipy import signal,interpolate

#from plotly.subplots import make_subplots
#import plotly.graph_objects as go
#import plotly.express as px

"""Adding Asthetics to the plot"""
plt.rc('font', size=15)          # controls default text sizes
plt.rc('axes', titlesize=25)     # fontsize of the axes title
plt.rc('axes', labelsize=15)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
plt.rc('ytick', labelsize=15)    # fontsize of the tick labels
plt.rc('xtick.major', size=14)    # size of the tick markers
plt.rc('ytick.major', size=14)    # size of the tick markers
plt.rc('xtick.minor', size=10)    # size of the tick markers
plt.rc('ytick.minor', size=10)    # size of the tick markers
plt.rc('legend', fontsize=15)    # legend fontsize
plt.rcParams['figure.figsize'] = [24, 12]  # set plotsize

global data_path, raw_db_filename, master_db_filename

"""Change these path and fileneames to your  dataset"""
data_path = r'/home/baparker/GitHub/Research/CCDs/Sophia_QE/20240607/'
# raw db has the dark measurements saved with 0 wavelength. Need to fix this in the future to avoid divide by zero error 
# raw_db_filename=r'D:\Nuvu_data\w18d10\qe_07112023\qe_07112023\MC_Charachterisation\May PTC Data\raw_MC_characterization_052023.csv'
# master db the dark corrected lamp data is coverted to photon count rate. Dark data is not recorded in this data frame to avoid divide by zero error. Need to fix this in the future and set a "wavelenght" for dark data for corresponding wavelength around which dark was taken.  
master_db_filename=r'/home/baparker/GitHub/Research/CCDs/Sophia_QE/master_MC_characterization_20240626.csv'

photodiode_response_path = r'/home/baparker/GitHub/Research/CCDs/Sophia_QE/photodiode_response.csv'
# BP Path that describes the photodiode A/W response at a given wavelength.

def get_photodiode_AW_response(wl):
    """
    Calculate the photodiode's A/W response at a given wavelength using interpolation.

    Parameters:
    wl (float): Input wavelength in nanometers.

    Returns:
    float: A/W response in Amperes per Watt.
    """
    # Read VUV Si response data from a CSV file
    VUV = pd.read_csv(photodiode_response_path)
    
    # Extract A/W response values from the data
    AW = pd.DataFrame(VUV, columns=['y'])
    AW['y'] = AW['y'].astype('float')
    
    # Extract wavelength values from the data
    lamb = pd.DataFrame(VUV, columns=['x'])
    lamb['x'] = lamb['x'].astype('float')
    
    # Define a set of quantile knots for spline interpolation
    knot_numbers = 20
    x_new = np.linspace(0, 1, knot_numbers + 2)[1:-1]
    q_knots = np.quantile(lamb, x_new) 
    
    # Perform spline interpolation using BSpline
    t, c, k = interpolate.splrep(lamb, AW, t=q_knots, s=1)
    photodiode_resfit = interpolate.BSpline(t, c, k)(wl)
    
    return photodiode_resfit

def currenttophotonrate(wl,current):
    """This function converts the measured signed photocurrent (in Amperes) of the Mcpherosn VUV photodiode measured using a picoammeter at a given wavelength (in nanometers) to photon count rate (photons/s) 

    Args:
        wl (float): wavelength of measured photocurrent in nanometers (nm)
        current (float): Signed Photocurrent measured from the VUV photodidoe in Amperes (A)

    Returns:
        photonrate (Astropy Quantity): returns photon rate in (photons/s). The ouput is an Astropy Quantity object with units of 1/s
    """
    
    """Import astropy constant and units package"""
    import astropy.constants as const
    from astropy import units as u
    """Remove the sign from photocurrent ad convert it to Astropy quanity object with units of Ampere"""
    picoAmpcurrent=-1*current*u.A
    """Get the photodiode A/W factor from measuremetns and assign it units of A/W"""
    AbyW=get_photodiode_AW_response(wl)*u.A/u.W
    """Calculate the picoammeter power uisng A/W """
    picoa_power_watt=(picoAmpcurrent/AbyW)
    """Calculate energy of a photon at a given wavelength """
    photon_energy=(const.h*const.c)/(wl*10**-9*u.m)
    """Convert picoammter power to photon count rate"""
    photonrate=picoa_power_watt/(photon_energy)
    # print(photonrate)
    #print(type(photonrate))
    return photonrate.decompose().value

def get_meta_data_from_filename(full_filename):
    """
    Extracts and returns metadata information from a given filename.

    Parameters:
    full_filename (str): The full path or name of the file.

    Returns:
    dict: A dictionary containing the extracted metadata information.
        - 'wavelength': Wavelength value rounded to 2 decimal places.
        - 'data_type': Type of data ('Dark', 'Xe_lamp', 'D2_lamp').
        - 'Slitsize': Size of the slit in microns.
        - 'Filter': Filter used in the experiment.
        - 'Ch1': Port information for Channel 1.
        - 'Ch2': Port information for Channel 2.
    """

    filename = os.path.basename(full_filename)

    # Check if the filename contains "cal_on" or "cal_off" and assign port information accordingly
    if filename.find("cal_on") > 0: 
        Ch1_port = "Main"
        Ch2_port = "Side"
    else: 
        Ch1_port = "Side"
        Ch2_port = "Main"

    # Check if the filename contains "dark" to determine the data type
    if filename.find("dark") > 0: 
        flag = "Dark"
        
        # Extract the slitsize from the filename
        loc = filename.find("slit")
        loc1 = filename.find("micron")
        slitsize = float(f"{filename[loc + 5:loc1]}")
        
        # Extract the filter information from the filename
        loc = filename.find("Filter")
        Filter = f"{filename[loc + 7:loc + 8]}"
        
        # Check if the filename contains wavelength information
        if filename.find("wl") > 0:
            loc = filename.find("wl")
            loc1 = filename.find("nm")
            wavelength = float(f"{filename[loc + 3:loc1]}")
        else:
            wavelength = 0.0

        return {
            'wavelength': [np.round(wavelength, 2)],
            'data_type': [flag],
            'Slitsize': [slitsize],
            'Filter': [Filter],
            'Ch1': [Ch1_port],
            'Ch2': [Ch2_port]
        }
    # Check if the filename contains wavelength information
    elif filename.find("wl") > 0:
        # Determine the data type based on lamp type
        if filename.find("Xe") > 0:
            flag = f"Xe_lamp"
        elif filename.find("D2") > 0:
            flag = f"D2_lamp"
        
        # Extract the slitsize from the filename
        loc = filename.find("slit")
        loc1 = filename.find("micron")
        slitsize = float(f"{filename[loc + 5:loc1]}")
        
        # Extract the filter information from the filename
        loc = filename.find("Filter")
        Filter = f"{filename[loc + 7:loc + 8]}"
        
        # Extract the wavelength information from the filename
        loc = filename.find("wl")
        loc1 = filename.find("nm")
        wavelength = float(f"{filename[loc + 3:loc1]}")

        return {
            'wavelength': [np.round(wavelength, 2)],
            'data_type': [flag],
            'Slitsize': [slitsize],
            'Filter': [Filter],
            'Ch1': [Ch1_port],
            'Ch2': [Ch2_port]
        }

def create_lamp_database(path):
    """
    Creates a pandas DataFrame containing data and metadata extracted from CSV files in the specified path.

    Parameters:
    path (str): The path where CSV files are located.

    Returns:
    pd.DataFrame: A DataFrame containing extracted data and metadata.
    """
    # Initialize an empty DataFrame with predefined columns
    df = pd.DataFrame(data=None, columns=['wavelength', 'data_type', 'Slitsize', 'Filter', 'Ch1', 'Ch2',
                                          'Ch1_total', 'Ch1_mean', 'Ch1_max', 'Ch1_min', 'Ch1_std', 'Ch1_median',
                                          'Ch2_total', 'Ch2_mean', 'Ch2_max', 'Ch2_min', 'Ch2_std', 'Ch2_median',
                                          'Filename', 'Expname'], dtype='float64')

    # Get a list of CSV files recursively in the specified path
    filelist = glob.glob(path + "*.csv", recursive=True)
    
    # Initialize variables to store mean dark values
    mean_dark_Ch1 = 0.0
    mean_dark_Ch2 = 0.0

    # Loop through each file in the list
    for idx, filename in enumerate(filelist):
        # Read data from the CSV file
        data = pd.read_csv(filename)

        # Check if the file is a dark calibration file
        if os.path.basename(filename).find("dark") > 0:
            mean_dark_Ch1 = np.mean(data.Ch1)
            mean_dark_Ch2 = np.mean(data.Ch2)
        else:
            # Subtract mean dark values from Ch1 and Ch2 data
            data.Ch1 = data.Ch1 - mean_dark_Ch1
            data.Ch2 = data.Ch2 - mean_dark_Ch2

        # Get metadata from the filename using the previously defined function
        #meta_data = pd.DataFrame(get_meta_data_from_filename(filename))

        # Append metadata to the DataFrame
        #df.loc[idx, meta_data.columns.to_list()] = meta_data.loc[0]

        # Compute statistics for Ch1 and Ch2 data
        data_stats = data.describe().transpose()

        # Update DataFrame with computed statistics
        df.loc[idx, ['Ch1_mean', 'Ch1_max', 'Ch1_min', 'Ch1_std']] = list(data_stats.loc[['Ch1'], ['mean', 'max', 'min', 'std']].values[0])
        df.loc[idx, ['Ch2_mean', 'Ch2_max', 'Ch2_min', 'Ch2_std']] = list(data_stats.loc[['Ch2'], ['mean', 'max', 'min', 'std']].values[0])
        df.loc[idx, ['Ch1_total', 'Ch1_median', 'Ch2_total', 'Ch2_median']] = [np.sum(data['Ch1']), np.median(data['Ch1']), np.sum(data['Ch2']), np.median(data['Ch2'])]

        # Store filename and experiment name
        df.loc[idx, ['Filename']] = os.path.basename(filename)
        df.loc[idx, ['Expname']] = os.path.basename(os.path.dirname(filename))
    df.to_csv(master_db_filename, index=None)
    print(f"Raw Database for the Monochromator charachterisation is save as: {master_db_filename}")
    return df

def update_lamp_database(master_df, master_db_filename=r'D:\Nuvu_data\w18d10\qe_07112023\qe_07112023\MC_Charachterisation\May PTC Data\master_df_validated_052023.csv'):
    """
    Updates a master DataFrame with additional calculated columns, saves it to a CSV file, and returns the updated DataFrame.

    Parameters:
    master_df (pd.DataFrame): The master DataFrame to be updated.
    master_df_path (str): The path to save the updated DataFrame as a CSV file.

    Returns:
    pd.DataFrame: The updated master DataFrame.
    """
    # Calculate various photon rates based on existing DataFrame columns. We drop the units from the current to photon rate conversion for convenience. 
    master_df['Ch1_mean_phrate'] = currenttophotonrate(master_df.wavelength.values, master_df.Ch1_mean.values) #photon rate in photons/second
    master_df['Ch1_median_phrate'] = currenttophotonrate(master_df.wavelength.values, master_df.Ch1_median.values)
    master_df['Ch1_min_phrate'] = currenttophotonrate(master_df.wavelength.values, master_df.Ch1_min.values)
    master_df['Ch1_max_phrate'] = currenttophotonrate(master_df.wavelength.values, master_df.Ch1_max.values)
    master_df['Ch1_std_phrate'] = np.abs(currenttophotonrate(master_df.wavelength.values, master_df.Ch1_std.values))
    master_df['Ch2_mean_phrate'] = currenttophotonrate(master_df.wavelength.values, master_df.Ch2_mean.values)
    master_df['Ch2_median_phrate'] = currenttophotonrate(master_df.wavelength.values, master_df.Ch2_median.values)
    master_df['Ch2_min_phrate'] = currenttophotonrate(master_df.wavelength.values, master_df.Ch2_min.values)
    master_df['Ch2_max_phrate'] = currenttophotonrate(master_df.wavelength.values, master_df.Ch2_max.values)
    master_df['Ch2_std_phrate'] = np.abs(currenttophotonrate(master_df.wavelength.values, master_df.Ch2_std.values))
    master_df['Ch2_std_phrate'] = currenttophotonrate(master_df.wavelength.values, master_df.Ch2_std.values)
    
    # Calculate ratio of median photon rates and their uncertainties
    if np.mean(master_df['Ch1_median_phrate']) > np.mean(master_df['Ch2_median_phrate']):
        master_df['RatioMbyS'] = master_df['Ch1_median_phrate'] / master_df['Ch2_median_phrate']
    else:
        master_df['RatioMbyS'] = master_df['Ch2_median_phrate'] / master_df['Ch1_median_phrate']
    
    ratio = []
    ratio_mean_error = []
   
    # Loop through rows and calculate ratio and its error
    for index, row in master_df.iterrows():
        if row['Ch1_median_phrate'] > row['Ch2_median_phrate']:
            try:
                ratio.append(row['Ch1_median_phrate'] / row['Ch2_median_phrate'])
                ratio_mean_error.append((row['Ch1_median_phrate'] / row['Ch2_median_phrate']) * np.sqrt(((row['Ch1_std_phrate'] / row['Ch1_median_phrate']) ** 2) + ((row['Ch2_std_phrate'] / row['Ch2_mean_phrate']) ** 2)))
            except ZeroDivisionError:
                ratio.append(0)
                ratio_mean_error.append(0)
        else:
            try:
                ratio.append(row['Ch2_median_phrate'] / row['Ch1_median_phrate'])
                ratio_mean_error.append((row['Ch2_mean_phrate'] / row['Ch1_median_phrate']) * np.sqrt(((row['Ch1_std_phrate'] / row['Ch1_median_phrate']) ** 2) + ((row['Ch2_std_phrate'] / row['Ch2_median_phrate']) ** 2)))
            except ZeroDivisionError:
                ratio.append(0)
                ratio_mean_error.append(0)

    master_df['RatioMbyS'] = ratio
    master_df['RatioMbyS_err'] = ratio_mean_error
    
    # Save the updated DataFrame to a CSV file
    master_df.to_csv(master_db_filename, index=None)
    print(f"Master Database for the Monochromator charachterisation is save as: {master_db_filename}")
    return master_df

def main():
    print("Reading raw data from the folder ")
    df=create_lamp_database(data_path)
    df=df[df.wavelength>0.0] # to remove dark data that has wavelength label of 0 nm
    print("Updating the raw data with photon rate calculations and caculating the ratio between Main and Sideport flux")
    master_df=update_lamp_database(df,master_db_filename) #already write the lamp database into a csv with filename in global varialble master_db_filename

if __name__ == "__main__":
    main()