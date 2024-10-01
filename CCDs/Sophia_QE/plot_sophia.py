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

dates = ['20240718','20240719','20240722','20240724', '20240726', '20240906']

#Good, Good, Water, Water warm and cold, Good, warm and test with lamp

wls = [245, 290, 295, 300, 300, 300]

in_folder = [['QE'], ['QE'], ['QE', 'warm'], ['QE', 'First'], ['QE'], ['Test']]

vmins =  [[0], [-100], [0,   -50],  [-500,  6500], [-1000],  [-300]]
vmaxes = [[8000],     [2000], [2250, 200], [27500, 8000], [15000], [150]]

temps = [['-85'], ['-85'], ['-85', '15'], ['20', '10'], ['-85'], ['5']]

for i, date in enumerate(dates):
    wl = wls[i]
    for j, folder in enumerate(in_folder[i]):
        base_file = base_dir + '/' + date + '/' + folder
        
        science_file = base_file + r'/science_{:.0f}_2.fits'.format(wl)
        bias_file = base_file + r'/bias_{:.0f}.fits'.format(wl)
        dark_file = base_file + r'/dark_{:.0f}.fits'.format(wl)
        
        science = fits.getdata(science_file)[0].astype('float32')
        dark = fits.getdata(dark_file)[0].astype('float32')
        bias = fits.getdata(bias_file)[0].astype('float32')
        
        if date == '20240718':
            science = dark
            dark = 0
        
        reduced = (science - bias) - (dark - bias)
        
        temp = temps[i][j]
                
        fig, ax = plt.subplots(layout='tight', figsize=(5,4))

        image = ax.imshow(reduced, vmin=vmins[i][j], vmax=vmaxes[i][j], cmap='inferno')#norm=LogNorm(vmin=vmins[i][j], vmax=vmaxes[i][j] ))
        
        ax.set_xlabel('Pixel')
        ax.set_ylabel('Pixel')
        ax.set_title('{}: {} nm'.format(date,wl))
        ax.tick_params(axis='both', direction='in', which='both')
        plt.colorbar(image, label='Counts', ax=ax, fraction=0.046, pad=0.04)
        
        patches = [mlines.Line2D([0], [0], marker='.', markersize=0, label=temp + '$^{{\circ}}$ C', ls='')]
        plt.legend(handles=patches, loc='upper left', handletextpad=-2)                
        fig.tight_layout()
        plt.savefig(save_dir + r'/Science_{}_{:.0f}_{}.png'.format(date, wl, temp),dpi = 500)
        plt.show()
        
        if date=='20240722':
            blur = 100 # pixels
            psf = Gaussian2DKernel(blur)

            convolved_image = convolve_fft(reduced, psf, boundary='fill')
            
            fig, ax = plt.subplots(layout='tight', figsize=(5,4))

            image = ax.imshow(convolved_image, vmin=vmins[i][j], vmax=vmaxes[i][j]*0.85, cmap='inferno')#norm=LogNorm(vmin=vmins[i][j], vmax=vmaxes[i][j] ))
            
            ax.set_xlabel('Pixel')
            ax.set_ylabel('Pixel')
            ax.set_title('Blurred {}: {} nm'.format(date,wl))
            ax.tick_params(axis='both', direction='in', which='both')
            plt.colorbar(image, label='Counts', ax=ax, fraction=0.046, pad=0.04)
            
            patches = [mlines.Line2D([0], [0], marker='.', markersize=0, label=temp + '$^{{\circ}}$ C', ls='')]
            plt.legend(handles=patches, loc='upper left', handletextpad=-2)                
            fig.tight_layout()
            plt.savefig(save_dir + r'/Science_Blurred_{}_{:.0f}_{}.png'.format(date, wl, temp),dpi = 500)
            plt.show()
                        






