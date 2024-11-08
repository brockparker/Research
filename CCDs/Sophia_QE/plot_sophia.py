#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 10:49:46 2024

@author: baparker
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft
import glob
import cv2
# BP Neccessary imports.

plt.rc('axes', labelsize=14)
plt.rc('figure', titlesize=30)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('legend', fontsize=12)
# BP Stylization parameters.

base_dir = r'/home/baparker/GitHub/Research/CCDs/Sophia_QE/Data'
save_dir = r'/home/baparker/GitHub/Research/CCDs/Sophia_QE/Photos'
#base_dir = r'C:\Users\Brock\Documents\Git\Research\CCDs\Sophia_QE\Data'
#save_dir = r'C:\Users\Brock\Documents\Git\Research\CCDs\Sophia_QE\Photos'
# BP Directories with data and to save to.

############## Plot specific dates to diagnose water issue ##############

dates = ['20240718','20240719','20240722','20240724', '20240726', '20240906']
# BP Dates to analyze, must be names of folders.
#Good, Good, Water, Water warm and cold, Good, warm and test with lamp

wls = [245, 290, 295, 300, 300, 300]
# BP Wavelenghts in nanometers, must be in image names, i.e. science_250.fits, bias_615.fits.

in_folder = [['QE'], ['QE'], ['QE', 'warm'], ['QE', 'First'], ['QE'], ['Test']]
# BP Folders inside of the main data folder (date).

vmins =  [[0], [-100], [0,   -50],  [-500,  6500], [-1000],  [-300]]
vmaxes = [[8000],     [2000], [2250, 200], [27500, 8000], [15000], [150]]
# BP Scales for linear plotting.

#temps = [['-85'], ['-85'], ['-85', '15'], ['20', '10'], ['-85'], ['5']]
# BP Manually set temperatures.

for i, date in enumerate(dates):
    # BP Loop over all dates in list.
    wl = wls[i]
    for j, folder in enumerate(in_folder[i]):
        # BP Some dates have multiple folders inside, loop over these.
        base_file = base_dir + '/' + date + '/' + folder
        # BP Create base file directory from information above.
        
        science_file = base_file + r'/science_{:.0f}_2.fits'.format(wl)
        bias_file = base_file + r'/bias_{:.0f}.fits'.format(wl)
        dark_file = base_file + r'/dark_{:.0f}.fits'.format(wl)
        # BP Define location for each file.
        
        science = fits.getdata(science_file)[0].astype('float32')
        dark = fits.getdata(dark_file)[0].astype('float32')
        bias = fits.getdata(bias_file)[0].astype('float32')
        # BP Read in each file without opening.
        
        temp = fits.open(science_file)[0].header['PI Camera Sensor Temperature Reading']
        # BP Get temperature from science file header.

        if date == '20240718':
            science = dark
            dark = 0
        # BP Issues with 20240718 darks, just set the dark to 0 instead.
        
        reduced = (science - bias) - (dark - bias)
        # BP Reduce the data for each individual wavelength.
                
        fig, ax = plt.subplots(layout='tight', figsize=(5,4))

        image = ax.imshow(reduced, vmin=vmins[i][j], vmax=vmaxes[i][j], cmap='inferno')#norm=LogNorm(vmin=vmins[i][j], vmax=vmaxes[i][j] ))
        
        ax.set_xlabel('Pixel')
        ax.set_ylabel('Pixel')
        ax.set_title('{}: {} nm'.format(date,wl))
        ax.tick_params(axis='both', direction='in', which='both')
        plt.colorbar(image, label='Counts', ax=ax, fraction=0.046, pad=0.04)
        
        patches = [mlines.Line2D([0], [0], marker='.', markersize=0, label=temp + '$^{{\circ}}$ C', ls='')]
        plt.legend(handles=patches, loc='upper left', handletextpad=-2)   
        # BP Create custom legend with temperature.             
        fig.tight_layout()
        plt.savefig(save_dir + r'/Science_{}_{:.0f}_{}.png'.format(date, wl, temp),dpi = 500)
        plt.show()
        # BP Plot data for each date.
        
        if date=='20240722':
            # BP Blur the data for 20240722 to try to recreate data from 20241719, before water incident.
            blur = 100 # pixels
            psf = Gaussian2DKernel(blur)
            convolved_image = convolve_fft(reduced, psf, boundary='fill')
            # BP Blur image
            
            fig, ax = plt.subplots(layout='tight', figsize=(5,4))

            image = ax.imshow(convolved_image, vmin=vmins[i][j], vmax=vmaxes[i][j]*0.85, cmap='inferno')#norm=LogNorm(vmin=vmins[i][j], vmax=vmaxes[i][j] ))
            
            ax.set_xlabel('Pixel')
            ax.set_ylabel('Pixel')
            ax.set_title('Blurred {}: {} nm'.format(date,wl))
            ax.tick_params(axis='both', direction='in', which='both')
            plt.colorbar(image, label='Counts', ax=ax, fraction=0.046, pad=0.04)
            
            patches = [mlines.Line2D([0], [0], marker='.', markersize=0, label=temp + '$^{{\circ}}$ C', ls='')]
            plt.legend(handles=patches, loc='upper left', handletextpad=-2)     
            # BP Create custom legend with temperature.             
            fig.tight_layout()
            plt.savefig(save_dir + r'/Science_Blurred_{}_{:.0f}_{}.png'.format(date, wl, temp),dpi = 500)
            plt.show()
            # BP Plot image.   

############## Plot good QE data and create movie ##############

# !!!!!!!!!!!!!! TODO Should update this once master biases/darks are created and data is compiled

save_dir = r'/home/baparker/GitHub/Research/CCDs/Sophia_QE/Photos/Video'

dates = ['20240927', '20241009', '20241010', '20241011', '20241024']
# BP Dates to analyze, must be names of folders.

folder = 'QE'
# Folder insdie of dates that contains raw data.

for i, date in enumerate(dates):
    # BP Loop over all dates in list.
    base_file = base_dir + '/' + date + '/' + folder
    # BP Create base file directory from information above.
        
    science_files = glob.glob(base_file + '/science_???_2.fits')
    bias_files = glob.glob(base_file + '/bias_???.fits')
    dark_files = glob.glob(base_file + '/dark_???.fits')
    
    science_files.sort(), bias_files.sort(), dark_files.sort() 
    # BP Find all data files and sort by wavelength. Verify all files correspond with correct wavlength and calibration files.

    for science_file, bias_file, dark_file in zip(science_files, bias_files, dark_files):
        # BP Loop over all files in each directory.
                
        science = fits.getdata(science_file)[0].astype('float16')
        dark = fits.getdata(dark_file)[0].astype('float16')
        bias = fits.getdata(bias_file)[0].astype('float16')
        # BP Read in each file without opening.
        
        wl = science_file[-10:-7]
        # BP Extract wavelength in nm from image name. May need to adjust if file naming scheme changes.
        
        temp = fits.getheader(science_file)['PI Camera Sensor Temperature Reading']
        # BP Get temperature from science file header.
        
        reduced = (science - bias) - (dark - bias)
        # BP Reduce the data for each individual wavelength.
                
        fig, ax = plt.subplots(layout='tight', figsize=(5,4))
    
        image = ax.imshow(reduced, cmap='inferno', vmin=0, vmax=40000)#norm=LogNorm(vmin=vmins[i][j], vmax=vmaxes[i][j] ))
        
        ax.set_xlabel('Pixel')
        ax.set_ylabel('Pixel')
        ax.set_title('{}: {} nm'.format(date,wl))
        ax.tick_params(axis='both', direction='in', which='both')
        plt.colorbar(image, label='Counts', ax=ax, fraction=0.046, pad=0.04)
        
        patches = [mlines.Line2D([0], [0], marker='.', markersize=0, label=temp + '$^{{\circ}}$ C', ls='')]
        plt.legend(handles=patches, loc='upper left', handletextpad=-2)   
        # BP Create custom legend with temperature.             
        fig.tight_layout()
        plt.savefig(save_dir + r'/Science_{}.png'.format(wl),dpi = 500) # !!!!!!! RESET TO 500 TODO
        plt.show()
        # BP Plot data for each date.
        
movie_images = sorted(glob.glob(save_dir + '/*.png'))
# BP Find all generated png images.

frame = cv2.imread(movie_images[0])
height, width, layers = frame.shape
# Get png properties.

video_name = save_dir + '/sophia_movie.avi'
video = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'avi1'), 10, (width,height))
# BP Create video object.

for image in movie_images:
    video.write(cv2.imread(image))
# BP Write all images to video object.

cv2.destroyAllWindows()
video.release()
# BP Save video and close video object.
