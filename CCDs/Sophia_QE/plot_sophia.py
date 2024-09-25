#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 10:49:46 2024

@author: baparker
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
from matplotlib.colors import LogNorm

plt.rc('axes', labelsize=14)
plt.rc('figure', titlesize=30)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('legend', fontsize=12)
# BP Stylization parameters.

# base_dir = '/home/GitHub/Research/CCDs/Sophia_QE/Data/'
base_dir = r'C:\Users\Brock\Documents\Git\Research\CCDs\Sophia_QE\Data'
save_dir = r'C:\Users\Brock\Documents\Git\Research\CCDs\Sophia_QE\Photos'

dates = ['20240718','20240719','20240722','20240724', '20240726', '20240906']

#Good, Good, Water, Water warm and cold, Good, warm and test with lamp

wls = [245, 290, 295, 300, 300, 300]

in_folder = [['QE'], ['QE'], ['QE', 'warm'], ['QE', 'First'], ['QE'], ['Test']]

vmins =  [[0],     [1250], [1500,  125], [20000,  6500], [15000], [100]]
vmaxes = [[-1000], [2000], [2000, 180],  [25000, 8000],  [17500], [150]]

for i, date in enumerate(dates):
    wl = wls[i]
    for j, folder in enumerate(in_folder[i]):
        base_file = base_dir + '\\' + date + '\\' + folder
        
        science_file = base_file + r'\science_{:.0f}_2.fits'.format(wl)
        bias_file = base_file + r'\bias_{:.0f}.fits'.format(wl)
        dark_file = base_file + r'\dark_{:.0f}.fits'.format(wl)
        
        science = fits.getdata(science_file)[0].astype('float32')
        dark = fits.getdata(dark_file)[0].astype('float32')
        bias = fits.getdata(bias_file)[0].astype('float32')
        
        reduced = (science - bias) - (dark - bias)
                
        fig, ax = plt.subplots(layout='tight', figsize=(5,6))

        image = ax.imshow(reduced, vmin=vmins[i][j], vmax=vmaxes[i][j], cmap='inferno')#norm=LogNorm(vmin=vmins[i][j], vmax=vmaxes[i][j] ))
        
        ax.set_xlabel('Pixel')
        ax.set_ylabel('Pixel')
        ax.set_title('{} {}: {} nm'.format(date,folder,wl))
        ax.tick_params(axis='both', direction='in', which='both')
        plt.colorbar(image, label='Counts', ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout(pad=1)
        plt.savefig(save_dir + r'test.png',dpi = 500)
        plt.show()
        
        #also gaussian blur the bad images to see if they look like the good images