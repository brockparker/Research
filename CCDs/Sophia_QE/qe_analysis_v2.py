"""
EMCCD_QQ_ANALYSIS - Python scripts and libararies for analysis of EMCCD QE measurement data taken with the UV-VIS Detector Characterisation Setup at Hamden UV/Vis Detector lab.
"""

__version__ = "0.0.2"

"""Importing libraries for OS commands and file systems"""
import os
import os.path
import glob

"""Importing libraries for numerical cacluations and data storage. """
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy import odr

"""Importing astropy modules for handling fits files, plotting and modeeling. """
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize,HistEqStretch,ZScaleInterval)
from astropy.modeling import models, fitting
from astropy.utils.exceptions import AstropyUserWarning

"""importing libraries for plotting."""
# from pylab import *
import matplotlib.pylab as plt
import matplotlib.colors as clr
import matplotlib.colors as mcolors
import random
import warnings

"""Updating Matplotlib options for plotting"""
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
plt.rcParams['figure.figsize'] = [18, 12]  # set plotsize
colors=['plum','green','gold','firebrick','dodgerblue','magenta','blue','lawngreen','cyan','navy','sienna','tomato','teal'] #default colors for plots

"""
Setting Global Variables to be used across various funcitons in this script
"""
global detector_name, conversion_gain
global dir_path, path_detector_metadata, path_coat_regions_metadata,path_analysis_regions_metadata,path_photodidode_database,dir_region_figs
global scan_name,coat_regions,analysis_regions,photodiode_data

""" 
Settings specific to detector being tested
"""
detector_name="w18d10"
conversion_gain=1.364 #e-/ADU

"""Location of the data files"""
#dir_path=r'D:\Nuvu_data\w18d10\qe_07112023\qe_07112023\qe5' ## folder with exposure, dark and bias taken during QE scans. For each wavelength we have one dark, one flat and one bias. 
dir_path=r'C:/Users/Brock/Documents/Git/Research/CCDs/Sophia_QE/20240607'
path_photodidode_database=r'D:\Nuvu_data\w18d10\qe_07112023\qe_07112023\qe5\Photodiode_db_qe_07112023_qe5.csv' ##filepath to database wiht measurements for photodiode at each wavelength. The current for the sideport photodidoe is logged through the entire exposure median dark and neduab count rate are extracted from the raw data using a seperate script. The photodidoe raw data is also located in the dir_path folder. 


"""Detector Meta data"""
path_detector_metadata=r"D:\Nuvu_data\w18d10\detector_metadata_w18d10\w18d10_metadata.csv" # specific information related to detector under test, eg detector name, lot number, number of coatings, processing history. 
path_coat_regions_metadata=r"D:\Nuvu_data\w18d10\detector_metadata_w18d10\w18d10_coat_regions_metadata.csv" # database of the different block coatings on the block coated detectors. 
path_analysis_regions_metadata=r"D:\Nuvu_data\w18d10\detector_metadata_w18d10\w18d10_analysis_regions_metadata.csv" # database of the different regions that will be extracted for QE analysis. These can be modified by modifying the csv file to set any arbitrary rectangular region for analyiss. In future version, we will have regions selected from pyds9 and astropy regions.  

"""File management: All figurees stored in figs folder in the same parent folder as the raw data. The regiosn folder is used to store images with the different regions that have been extracted for analysis"""
dir_region_figs=os.path.join(dir_path,r'figs\regions')
if os.path.exists(os.path.dirname(dir_region_figs))==False: #check if folder exists 
    os.mkdir(os.path.dirname(dir_region_figs)) # not create the folder 
if os.path.exists(dir_region_figs)==False: #check if folder exists 
    os.mkdir(dir_region_figs)  # not create the folder 

"""Name of the scan for creating assocaited data product names using the folder name in which the raw data is stored. These fodlernames are created delibrately to organize the datasets."""
scan_name=os.path.basename(dir_path)

# import all dataabases as pandas tables 
detector_metadata=pd.read_csv(path_detector_metadata)
coat_regions=pd.read_csv(path_coat_regions_metadata,skiprows=2) # skiprows to skip header rows) 
analysis_regions=pd.read_csv(path_analysis_regions_metadata,skiprows=2)# skiprows to skip header rows) 
photodiode_data=pd.read_csv(path_photodidode_database,index_col=0) # getting photodiode data

def generate_qe_db(dir_path):
    """This fuction crawls through the folder in which raw data from QE measurements is store and creates a database using the metadata in the fileanmes of the fits files, the log of the scan and the fits header. Theree dabases are creatd and stored as csv files for all images, bais images only, flat exposure images, and dark images.  
    it compares the filenames in the folder to the fileneames in the log to ensure that all scan files are present. 

    Args:
        dir_path (path): Path to the fodler wehre the raw data from QE scans is store.d 
        scan_log (csv file): Each scan folder with dir_path has a scan_log.csv file that identifies which image  
    Returns:
        images_db (pandas table): pandas datbase with information of all fits images (dark, bias, exposure)
        exposreu_db (pandas table): substed of the images_db created for convenience conatineing a list of flat exposure images with scan metadata, like associated wavelenght, filename etc. 
        dark_db (pandas table): substed of the images_db created for convenience conatineing a list of dark images with scan metadata, like associated wavelenght, filename etc. 
        bias_db (pandas table): subste of the images_db created for convenience conatinging only bias images with scan metadata, like assocaited wavelenght, filename etc.
        
    """    
    
    #get list of all files in the path folder
    dir_path_list=glob.glob(dir_path)

    print(f"Going through dir: {dir_path}")
    log_file=r'scan_log.csv' # name of the log file assocaited with each scan created by the monochromator scan script
    log_path=os.path.join(dir_path,log_file) #full path to the log file 
    images_db=pd.read_csv(log_path,index_col=0) # reading the log into a database 
    file_list=glob.glob(os.path.join(dir_path,'*.fits.Z')) #getting list of all files with the fits.Z extentions created by the nuvu acquisition scripts. 
    if len(file_list)!=np.shape(images_db)[0]: # check if number of fits files in the fodler are same as the number of fits files in the scan. 
        print(f"The number of files in folder is {len(file_list)} While the number of files in database are {np.shape(images_db)[0]}.")
        start_index= len(file_list)-np.shape(images_db)[0] #the file numbers may differ because of some test images taken before the scan was started. The start index is updated to start by skiping the first few files. 
        file_list=file_list[start_index:] # updating the file list using the start index. 
        images_db['filenames']=file_list # storing the filenlist in the database this stores the full path of the files in the database 
    else: 
        images_db['filenames']=file_list#  storing the filenlist in the database this stores the full path of the files in the database 

    # db_path=str(Path(dir_path).parents[0])
    db_path=dir_path # setting the databse torage path to same as the dir_path
    #extrating exposure images only, bias only, dark only images to a separate databases 
    exposure_db=images_db[images_db['imtype']=='Exposure'].reset_index() 
    bias_db=images_db[images_db['imtype']=='Bias'].reset_index()
    dark_db=images_db[images_db['imtype']=='Dark'].reset_index()    
    #saving the daatabses in the same  folder as raw data. 
    images_db.to_csv(os.path.join(db_path,'images_db.csv'))
    exposure_db.to_csv(os.path.join(db_path,'flat_db.csv'))
    bias_db.to_csv(os.path.join(db_path,'bias_db.csv'))
    dark_db.to_csv(os.path.join(db_path,'dark_db.csv'))
    return images_db,exposure_db,bias_db,dark_db #returnign the databases as pandaas tables. 

class extract_region:
    """Extract region creates region objects taht are used for storing metadata of regions extracted from the fits file for anlaysis. We use astropy cutouts for extracting regions in the fits files. the metadata tstored in the regions object is a helfpul tool when creating these cutouts. 
    """
    def __init__(self, xmin,xmax,ymin,ymax,name='img',color='green'):
        """The properties of each region object is stored here. Based on four concers fo the xmin, xmax, ymin, ymax fo the region to be extrated other properties like the center positions, shape, size etc are computed. Name and color atrributes are used to store 

        Args:
            xmin (float): x coordinate value of the bottom corner of the region. Origin pixel on the image as reference 
            xmax (float): x coordinate value of the top  corner of the region. Origin pixel on the image as reference 
            ymin (float): y coordinate value of the bottom corner of the region. Origin pixel on the image as reference
            ymax (float): y coordinate value of the top  corner of the region. Origin pixel on the image as reference 
            name (str, optional): name of the region for unique identification for e.g. unique name for a coaitng region or serial number an analysis region. Defaults to 'img'.
            color (str, optional): name of the region that can be used for plottingg. Use matplotlib color names. Defaults to 'green'.
        """
        self.xmin = int(xmin)
        self.xmax = int(xmax)
        self.ymin = int(ymin)
        self.ymax = int(ymax)
        self.centx=int((xmax+xmin)/2)
        self.centy=int((ymax+ymin)/2)
        self.sizex=int((xmax-xmin))
        self.sizey=int(ymax-ymin)
        self.centpos= (self.centx,self.centy)
        self.size=(self.sizey,self.sizex)
        self.color=color
        self.name=name
    # def describe():
    #     print(self.xmin)
        
def get_data(wl,exptime,coat,Region,center_position,size,Signal,Std,qe,image_fname='',bias_fname=''):
    # import datetime
    """Helper tool to format each row appended to a pandas dataframe. The function creates a pyton dictionary that can be appended to a pandas dataframe.  
    Args:
        wl (float): wavelength (nm)
        exptime (float): exposure time (seconds)
        coat (str): coating on the region 
        Region (str): region name
        center_position (tuple): tuple with x and y position of the geometeric center of the region (pixels)
        size (tuple): x and y dimensions of region (pixels)
        Signal (float): Counts in (ADU)
        Std (float): STD of the counts in (ADU)
        qe (float): QE (dimensionless ratio)
        image_fname (str, optional): full path of the flat exposure image. Defaults to ''.
        bias_fname (str, optional): full path of the bias exposure image used for correction. Defaults to ''.

    Returns:
        python dictionary: all metadata collected into  python discitin for injesion into pandas Dataframe 
    """
    data = {
        'wl': [],
        'exptime':[],
        'coat':[],
        'Region': [],
        'center_position': [],
        'size': [],
        'Signal':[],
        'Std':[],
        'QE':[],
        'image_fname':[],
        'bias_fname':[]
        }
    data['wl'].append(wl)
    data['exptime'].append(exptime)
    data['coat'].append(coat)
    data['Region'].append(Region)
    data['center_position'].append(center_position)
    data['size'].append(size)
    data['Signal'].append(Signal)
    data['Std'].append(Std)
    data['QE'].append(qe)
    data['image_fname'].append(image_fname)
    data['bias_fname'].append(bias_fname)
    return data

def get_reduced_image_wl(wl,bias_db,dark_db,Exposure_db): 
    biasfile=0
    print(f'Analysing data for {wl} nm')
    b_filename=bias_db[bias_db['wl']==wl].filenames.values[biasfile]
    imageb = np.array(fits.getdata(b_filename))
    ##modified from j to i index for one data set
    imageb=imageb*1.0
    d_filename=dark_db[dark_db['wl']==wl].filenames.values[0]   
    imaged= np.array(fits.getdata(d_filename))

    f_filename=Exposure_db[Exposure_db['wl']==wl].filenames.values[0]
    imagef = np.array(fits.getdata(f_filename)) 
    imagef=imagef##modified from j to i index for one data set
    imagef = imagef*1.0
    image= imagef-np.mean(imageb)-(imaged-np.mean(imageb))
    return image

def plot_regions_cutout(regions_list,wl_implot,bias_db,dark_db,Exposure_db): 
    global detector_name
    # wl_implot=170
    imagef=get_reduced_image_wl(wl_implot,bias_db,dark_db,Exposure_db)
    imagef = imagef*1.0
    fig,ax=plt.subplots(figsize=(60,20))
    interval = MinMaxInterval()
    vmin, vmax = interval.get_limits((imagef))
    # Create an ImageNormalize object using a SqrtStretch object
    # norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())
    #norm_max=vmax
    #norm_min=vmin
    
    norm_max=10000
    norm_min=1000
    norm=ImageNormalize(vmin=norm_min, 
                        vmax=norm_max, 
                        #stretch=SqrtStretch()
                        )
    im=ax.imshow((imagef),cmap='magma',norm=norm,origin='lower')
    cbar=fig.colorbar(im,extend='max')
    for idx,region in enumerate(regions_list): 
        position = region.centpos
        size = region.size    # pixels
        cutoutf1=Cutout2D(imagef, position, size,)
        cutoutf1.plot_on_original(color=colors[idx],label=region.name)
    ax.legend()
    ax.set_xlim(1082,2*1072)
    plt.grid(which='both',alpha=0.2)
    plt.title(f"{detector_name}@{wl_implot} nm, {scan_name}, exptime={Exposure_db.Exp_time[0]} seconds\n Analysis: Bias and Dark corrected")
    ax.set_xlabel("x pixels")
    ax.set_xlabel("y pixels")
    fig_filename=os.path.join(f"{dir_region_figs}",f"{scan_name}_{detector_name}_{wl_implot}nm.png")
    plt.savefig(fig_filename)

def run_qe_analysis(bias_db,dark_db,images_db,Exposure_db):
    """This scripts loops over ther region database and extracts the different regions for QE analysis form fits images. For each wavelength the regions are extracted and reduced by subracting the bias. For each region we caculate the mean flux and noise (std of the bias subracted signal). QE is caculated using 

    Args:
        bias_db (pandas table): database of all bias images for the scan 
        dark_db (pandas table): database of all dark images for the scan
        images_db (pandas table): database of all flat images for the scan
        Exposure_db (pandas table): database of all images for the scan
    """
    regions_list=[] #list to store a list of regions
    for index,row in analysis_regions.iterrows(): #iterate over rows of the pandas table 
        regions_list.append(extract_region(row.x1,row.x2,row.y1,row.y2,name=row.region_name,color=colors[index])) 
    # 
    # region_plot_wl=180
    # plot_regions_cutout(regions_list,region_plot_wl,bias_db,dark_db,Exposure_db)
    wllist=images_db['wl'].unique()

    biasfile=0 # which of the two bias files to use (pre or post)-- only one for the recent scans
    for idxb, wl in enumerate(wllist): #iterate over all wavelengths 
        #plots the full image with cutout location for each region and save it in a fodler. 
        plot_regions_cutout(regions_list,wl,bias_db,dark_db,Exposure_db)
        # read the bias and flat fits image data input numpy array
        print(f'Analysing data for {wl} nm')
        b_filename=bias_db[bias_db['wl']==wl].filenames.values[biasfile] #selecting bias file from the database
        imageb = np.array(fits.getdata(b_filename)) #numpy array
        imageb=imageb*1.0#enforcing all numbers to be float 
        # d_filename=dark_db[dark_db['wl']==wl].filenames.values[0]   #selecting dark file from the database
        # imaged= np.array(fits.getdata(d_filename))
        f_filename=Exposure_db[Exposure_db['wl']==wl].filenames.values[0] #selecting flat file from the database
        imagef = np.array(fits.getdata(f_filename)) #numpy array
        imagef = imagef*1.0 #enforcing all numbers to be float        
        reg=1 #setting region counter. 
        for region in regions_list: # iterating over the list of regions 
            regsize = (region.sizey,region.sizex) # size of extaction region as a tuple
            center_pos= region.centpos # center of extraction region as a tuple 
            # print(center_pos)
            imagef_sec=Cutout2D(imagef, center_pos, regsize).data # extract the region from bias image
            imageb_sec=Cutout2D(imageb, center_pos, regsize).data # extract the region form dark image
            # imaged_sec=Cutout2D(imaged, center_pos, regsize).data
            
            # print(np.mean(imagef_sec))
            # print(images_db.Exp_time[0])
            sig= np.std(imagef_sec-np.mean(imageb_sec)) #cacualte standard devaition of the signal in the regon after subrating mean bias 
            mf = np.mean(imagef_sec-np.mean(imageb_sec)) #cacualte mean signal in the region after subrating mean bias 
            
            # QE= (mean counts in the reion x convert signal form from ADu to e-/exposure time)/estiamted photon count rate
            # Exposure time is the same for all exposures  we take the Exptime from the first file in the images_db.
            qe=((mf*conversion_gain)/images_db.Exp_time[0])/(photodiode_data[photodiode_data.Wavelength==wl].Photodidode_data_counts_per_pix.values[0])
            #print(qe)
            # coat=whichcoat(center_pos,regsize)
            
            # storing the caculated values in the pandas dataframe 
            if idxb==0:
                # creates a new dataframe using the first entry  
                region_data =pd.DataFrame(get_data(wl=wl,exptime=images_db.Exp_time[0],coat="",Region=region.name,center_position=center_pos,size=regsize,Signal=mf, Std=sig, qe=qe ,image_fname=f_filename,bias_fname=b_filename))
            else:
                # appends row for a regionin the dataframe
                temp_data =pd.DataFrame(get_data(wl=wl,exptime=images_db.Exp_time[0],coat="",Region=region.name,center_position=center_pos,size=regsize,Signal=mf, Std=sig,qe=qe,image_fname=f_filename,bias_fname=b_filename))
                region_data=pd.concat([region_data,temp_data],ignore_index=True)
                
            reg=reg+1
    #saves the extracted data to file 
    region_data.to_csv(os.path.join(dir_path,f'regiondata_{scan_name}.csv')) 
    plot_data=1
    if plot_data==1: 
        #plotting the countrate and saving the plot to file 
        fig,ax=plt.subplots(figsize=(18,12))
        region_data['Region']=region_data['Region'].astype(str)
        for idx,regnum in enumerate(region_data.Region.unique()): #itereate over different unique region to plot for each region
            ax.scatter(region_data[region_data.Region==regnum].wl,region_data[region_data.Region==regnum].Signal/images_db.Exp_time[0],label=f'{regnum}',color=colors[idx],s=15)
        ax.legend(title="Regions")
        ax.set_xlabel("wavelength(nm)")
        ax.set_ylabel("count rate (e-/s/pix)")
        ax.set_yscale('log')
        ax.set_xlim(150,600)
        plt.title(f"Mean count rate (ph/s/pix) in different regions vs wavelength")
        # ax.set_xlim(160,420)
        # ax.set_ylim(0,2500)
        ax.grid(color='Grey',alpha=0.2,linestyle='-')    
        plt.savefig(os.path.join(dir_path,'countrate.png')) #plot saved in he same path as raw data (Chnage this to fig folder for future)
        print(f"Count rate vs wavelength plot saved as:{os.path.join(dir_path,'countrate.png')}")
        
        #plotting the QE vs wavelength and saving the plot to file 
        fig,ax=plt.subplots(figsize=(18,12))
        region_data['Region']=region_data['Region'].astype(str)
        for idx,regnum in enumerate(region_data.Region.unique()):
            ax.scatter(region_data[region_data.Region==regnum].wl,region_data[region_data.Region==regnum].QE,label=f'{regnum}',color=colors[idx],s=15)
        ax.legend(title="Regions")
        ax.set_xlabel("wavelength(nm)")
        ax.set_ylabel("QE")
        #ax.set_yscale('log')
        ax.set_xlim(150,600)
        plt.title(f"QE for different regions vs wavelength")
        # ax.set_xlim(160,420)
        ax.set_ylim(0,1)
        ax.grid(color='Grey',alpha=0.2,linestyle='-')    
        plt.savefig(os.path.join(dir_path,'qe.png')) #plot saved in he same path as raw data (Chnage this to fig folder for future)
        print(f"QE Plots saved as:{os.path.join(dir_path,'qe.png')}")

def main(): 
    # generate the databases for the fits files     
    images_db,Exposure_db,bias_db,dark_db=generate_qe_db(dir_path)
    # run the QE analysis 
    run_qe_analysis(bias_db,dark_db,images_db,Exposure_db)
    
    
if __name__ == "__main__":
    main()
    
    
""" Deprecated helper fulctions"""
# def set_regions(Region,coat,region_centx,region_centy,region_wx,region_wy):
#     """_summary_

#     Args:
#         Region (_type_): _description_
#         coat (_type_): _description_
#         region_centx (_type_): _description_
#         region_centy (_type_): _description_
#         region_wx (_type_): _description_
#         region_wy (_type_): _description_

#     Returns:
#         _type_: _description_
#     """
#     {'Region','Coat','region_centx','region_centy','region_wx','region_wy'}
#     data = {
#         'Region': [],
#         'coat':[],
#         'region_centx': [],
#         'region_centy': [],
#         'region_wx': [],
#         'region_wy':[]
#         }
#     data['Region'].append(Region)
#     data['Coat'].append(coat)
#     data['region_centx'].append(region_centx)
#     data['region_centy'].append(region_centy)
#     data['region_wx'].append(region_wx)
#     data['region_wy'].append(region_wy)
#     return data

# def get_data(wl,exptime,coat,Region,center_position,size,Signal,Std,qe,image_fname='',bias_fname=''):
#     import datetime
#     data = {
#         'wl': [],
#         'exptime':[],
#         'coat':[],
#         'Region': [],
#         'center_position': [],
#         'size': [],
#         'Signal':[],
#         'Std':[],
#         'QE':[],
#         'image_fname':[],
#         'bias_fname':[]
#         }
#     data['wl'].append(wl)
#     data['exptime'].append(exptime)
#     data['coat'].append(coat)
#     data['Region'].append(Region)
#     data['center_position'].append(center_position)
#     data['size'].append(size)
#     data['Signal'].append(Signal)
#     data['Std'].append(Std)
#     data['QE'].append(qe)
#     data['image_fname'].append(image_fname)
#     data['bias_fname'].append(bias_fname)
#     return data