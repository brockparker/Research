"""
QE_Phodiode_analysis - Python script for analysis of photocurrent measurements taken along with QE chaaracterisation exposures with UV-VIS Detector Characterisation Setup at Hamden UV/Vis Detector lab. This scripts generates a database that is used for deriving the QE from EMCCD image data at different wavelengths.
"""

__version__ = "0.0.2"


import numpy as np
import os
import glob
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import astropy.constants as const
from astropy import units as u
stop

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
plt.rcParams['figure.figsize'] = [18, 22]  # set plotsize

""" MC photiodiode scan holds database with the data from the charachterisation of the monochromator"""
global mc_photodiode_scan_path,qe_photodiode_scan_path,photo_diode_AbyW_response_path
mc_photodiode_scan_path=r'D:\Nuvu_data\w18d10\qe_07112023\qe_07112023\MC_Charachterisation\May PTC Data\master_df_validated_052023.csv'
qe_photodiode_scan_path=r'D:\Nuvu_data\w18d10\qe_07112023\qe_07112023\qe5'
photo_diode_AbyW_response_path=r'D:\Monochromator_data\vuvSiresponse_UV_PhotonResponse_plot-data.csv'

def get_photodide_AW_response(wl,path=photo_diode_AbyW_response_path): 
    """ Get the Amperes per Watts conversion factor for the photodidoe for the given waelenght from a csv file
    Args:
        wl (float): Input wavelength in nm
        path (str): Path to the file that contains the AbyW charachteristics of the Photodiode
    Returns:
        float: Ouput A/W in Amperes per Watt for the input wavelenght
    """

    VUV = pd.read_csv(path)
    AW = pd.DataFrame(VUV, columns=['y'])
    AW['y'] = AW['y'].astype('float')
    
    lamb = pd.DataFrame(VUV, columns=['x'])
    lamb['x'] = lamb['x'].astype('float')
    # if (wl>max(lamb.x))|(wl<min(lamb.x)): 
    #     print("Input wavlength greater than Max available data. The output is not be reliable fit")
    #     return 0        
    #spline
    knot_numbers = 20
    x_new = np.linspace(0, 1, knot_numbers+2)[1:-1]
    q_knots = np.quantile(lamb, x_new) 
    t,c,k = interpolate.splrep(lamb, AW, t=q_knots, s=1)
    photodiode_resfit = interpolate.BSpline(t,c,k)(wl)
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
    
    """Remove the sign from photocurrent ad convert it to Astropy quanity object with units of Ampere"""
    picoAmpcurrent=-1*current*u.A
    """Get the photodiode A/W factor from measuremetns and assign it units of A/W"""
    AbyW=get_photodide_AW_response(wl)*u.A/u.W
    """Calculate the picoammeter power uisng A/W """
    picoa_power_watt=(picoAmpcurrent/AbyW)
    """Calculate energy of a photon at a given wavelength """
    photon_energy=(const.h*const.c)/(wl*10**-9*u.m)
    """Convert picoammter power to photon count rate"""
    photonrate=picoa_power_watt/(photon_energy)
    return photonrate

def generate_MbyS_lookup(mc_photodiode_scan=r'D:\Nuvu_data\w18d10\qe_07112023\qe_07112023\MC_Charachterisation\May PTC Data\master_df_validated_052023.csv'): 
    # """ The Monochrmator charachtersiation database consists of data collected during charachtersiation of the monochromator with photodiodes at both port. This script #extacts the wavelength vs Mainport/Sideport ratio from this datbase into a lookup table (a Pandas dataframe). The MbyS ratio is unitless ratio of measured #photocurents from the photodiode at the two prts.  

    #     Args:
    #         mc_photodiode_scan (regexp, optional): Path to the database with measured charachtersistics of the monochromator output (eg. Mean current, median current etc. vs wavelength). Defaults to r'D:\Nuvu_data\w18d10\qe_07112023\qe_07112023\MC_Charachterisation\May PTC Data\master_df_validated_052023.csv'.
    #     #Returns:
    #         padndas dataframe: pandas dataframe with two colummns wavelength and the MbyS value for the corresponding wavelenght """
    master_df=pd.read_csv(mc_photodiode_scan_path) # read database 
    master_df=master_df[master_df.data_type=='D2_lamp'] # extract data taken with D2 lamp
    temp_table=master_df[['wavelength', 'RatioMbyS']] # get the wavelength vs MbySideRatio 
    temp_table=master_df[master_df.wavelength>0] # filter out wavelenght less than 0 (0 used for dark measurements)
    lookup_table=pd.DataFrame() #setupa lookup table dataframe
    lookup_table=temp_table.groupby('wavelength')['RatioMbyS'].median() #the database has entries from multiple scans take median value for each wavelenght and create a signle valued lookup table for wavelenght vs RatioMbyS. 
    lookup_table=lookup_table.reset_index() # Reinded and return the look_up table
    return lookup_table

def get_qe_meta_data_from_filename(full_filename):
    """Reads the filename of the photodiode data file generated from the QE scan code. The files names have metadata infomration including the filter type, image type, wavelength etc. This fucntion parses the filenames and retuns the meta from the filname as a python dictionary. 

    Args:
        full_filename (str): filename from the scan which has to be paresed. 

    Returns:
        dictionary: dictionary contain the metadata generated from the input filename
    """
    filename=os.path.basename(full_filename)
    if filename.find("Dark")>0: 
        data_type="Dark"
    else: 
        data_type="Image"
    
    if filename.find("f1")>0: 
        Filter="f1"
    elif filename.find("f2")>0: 
        Filter="f2"
    elif filename.find("f3")>0: 
        Filter="f3"
    elif filename.find("f4")>0: 
        Filter="f4"
    elif filename.find("f2")>0: 
        Filter="f5"
    else: 
        print("No filter information in filename. Check the file naming and file list.")
        return()

    loc=filename.find("nm")
    wavelength= float(f"{filename[loc-3:loc]}")
    
    imno=int(f"{filename[loc+3:-4]}")
    filename=os.path.basename(filename)
    return {'wavelength':[np.round(wavelength,2)],'data_type':[data_type],'Filter':[Filter],'imno':[imno],'filename':[filename],'fullfilename':[full_filename]}


def generate_qe_photodiode_db(mc_photodiode_scan_path,qe_photodiode_scan_path):
    """This function takes the paths for the Monochromator charachterisation dataset and the path for the photodidoe data from the QE scan. Based on the ratio between the flux between Sideport and Mainport taken during the monochromator charachtersiation, the script scales the measurement done at ths dideport during QE scans to predict the flux on the detector. The output is stored as a lookup table. This script has to be run for each QE measurement only once to generate the lookup table that is saved as a csv file. 

    Args:
        mc_photodiode_scan_path (str): _description_
        qe_photodiode_scan_path (str): _description_
    """
    from astropy import units as u 
    from astropy import constants as const
    
    MbyS_lookup=generate_MbyS_lookup(mc_photodiode_scan_path)
    scan_name=f'{os.path.basename(os.path.dirname(qe_photodiode_scan_path))}_{os.path.basename(qe_photodiode_scan_path)}'
    qe_df=pd.DataFrame(data=None,columns=['wavelength','data_type','Filter','imno','filename','fullfilename'],dtype='float64')
    dir_path_list=glob.glob(qe_photodiode_scan_path)
    filename='ptc_log.csv'
    
    print(f"Going through dir: {qe_photodiode_scan_path}")
    file_list=glob.glob(os.path.join(qe_photodiode_scan_path,'picoa*.csv'))
    #print(file_list)
    # print(file_list)
    for idx,filename in enumerate(file_list): 
        pd_data =pd.read_csv(filename)
        meta_data=pd.DataFrame(get_qe_meta_data_from_filename(filename))
        #append metadata to pandas dataframe
        qe_df.loc[idx,meta_data.columns.to_list()]=meta_data.loc[0]
    print(qe_df)
    
    
    wl_list=qe_df.wavelength.unique()
    Filter='f1'
    qe_df_f1=qe_df[qe_df['Filter']==Filter]

    median_dark=[]
    median_data=[]
    pickoff_data=[]
    mainport_estimate=[]
    pickoff_data_per_pix=[]
    mainport_estimate_per_pix=[]
    pickoff_area= 11*u.mm*u.mm #in mm^2 (from sideport anaylysis)
    main_photo_area=100*u.mm*u.mm # in mm^2 photodiode 10x10 mm filled completely. 
    # MbySarea=main_photo_area/pickoff_area #photodiode pickoff area =11 sq mm, #photodiode area =100 sq mm
    EMCCD_pixel_area=(13.0*13.0*u.micron*u.micron).decompose()#in mm^2/px (each pixel is 13.0 microns)
    npixels_mainport=((main_photo_area)/EMCCD_pixel_area).decompose()
    npixels_sideport= ((pickoff_area)/EMCCD_pixel_area).decompose()
    
    for idx,wl in enumerate(wl_list): 
        MbyS_wl_ratio=MbyS_lookup[MbyS_lookup.wavelength==wl].RatioMbyS.values[0]
        print(MbyS_wl_ratio)
        file_names=qe_df_f1[qe_df_f1.wavelength==wl].fullfilename.values
        for file_name in file_names: 
            if file_name.find("Dark")>0: 
                dark_data =pd.read_csv(file_name)
                dark_cal=np.median(dark_data.Ch1)
                median_dark.append(dark_cal)
                print(f'median dark for {wl}nm={currenttophotonrate(wl,median_dark[idx]).decompose().value}')
            elif file_name.find("Exposure")>0: 
                image_data=pd.read_csv(file_name)
                median_cal=np.median(image_data.Ch1)
                median_data.append(median_cal)
                # print(median_cal)
                pickoff_data.append(currenttophotonrate(wl,median_cal-dark_cal).decompose().value)
                mainport_estimate.append(MbyS_wl_ratio*(currenttophotonrate(wl,median_cal-dark_cal).decompose().value))
                pickoff_data_per_pix.append((((1/npixels_sideport)*currenttophotonrate(wl,median_cal-dark_cal))).decompose().value)
                mainport_estimate_per_pix.append((MbyS_wl_ratio*((1/npixels_mainport)*currenttophotonrate(wl,median_cal-dark_cal))).decompose().value)
                print(f'median flux for {wl}nm={currenttophotonrate(wl,median_cal-dark_cal).decompose().value}')

    # emccd_pixels=pickoff_area/EMCCD_pixel_area
    fig,ax=plt.subplots(figsize=(18,12))
    plt.scatter(wl_list,mainport_estimate_per_pix,color='Green',marker='x',label='mainport')  
    plt.scatter(wl_list,pickoff_data_per_pix,color='Blue',marker='x',label='sideport') 
    # ax.legend(title="Regions")
    ax.set_xlabel("wavelength(nm)")
    ax.set_ylabel("Photons/s/pix")
    ax.set_yscale('log')
    ax.set_xlim(150,600)
    # ax.set_ylim(1E+3,2E+5)
    plt.title(f"Photodiode median count rate/pix estimated at detector \n scan name= {scan_name}")
    # ax.set_xlim(160,420)
    # ax.set_ylim(0,2500)
    ax.grid(color='Grey',alpha=0.2,linestyle='-')

    plt.legend()

    pd_qe_scan_database=pd.DataFrame({'Wavelength':[],'Filter':'','imno':[],'Photodidode_data_counts_per_pix':[],'Mainport_est_counts_per_pix':[],"Expname":'','Sideport_Filename':'','MbySRatio_Filename':''})
    pd_qe_scan_database['Wavelength']=wl_list
    pd_qe_scan_database['Expname']=scan_name
    pd_qe_scan_database['Sideport_Filename']=qe_df_f1[(qe_df_f1.data_type=='Image')].fullfilename.values
    pd_qe_scan_database['imno']=qe_df_f1[(qe_df_f1.data_type=='Image')].imno.values
    pd_qe_scan_database['MbySRatio_Filename']=mc_photodiode_scan_path
    pd_qe_scan_database['Photodidode_data_counts_per_pix']=pickoff_data_per_pix
    pd_qe_scan_database['Mainport_est_counts_per_pix']=mainport_estimate_per_pix
    pd_qe_scan_database['Filter']=Filter
    db_path=os.path.join(qe_photodiode_scan_path,f'Photodiode_db_{scan_name}.csv')
    pd_qe_scan_database.to_csv(db_path)
    print(db_path)


def main():
    """Run the funciton below in the main string. Modifiy the path names in the global parameters at the beginning of this cript.
    """
    generate_qe_photodiode_db(mc_photodiode_scan_path,qe_photodiode_scan_path)
    
if __name__ == "__main__":
    main()