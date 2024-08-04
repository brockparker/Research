# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 14:10:10 2024

@author: Brock
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

base_dir = r'C:\Users\Brock\Documents\Git\Research\CCDs\Fermilab\LTA_Data\\'

file = base_dir + 'image_lta_07.fz'
new_file = base_dir + 'test_image_lta_07.fits'

hdul = fits.open(file)

def function(x):
    return x[0]**2 + x[1]

for amplifier in range(4):
    
    data = hdul[amplifier + 1].data
    hdr = hdul[amplifier + 1].header
    
    ncol = int(hdr['NCOL'])
    ccdncol = int(hdr['CCDNCOL'])
    ccdnpres = int(hdr['CCDNPRES']) + 1
    
    nover = int(ncol - ccdncol/2 - ccdnpres)
    
    overscan = data[:, -nover:]
    
    plt.imshow(data)
    plt.show()
        
# how to stitch together??? which quadrant goes where
    
hdul.writeto(new_file, overwrite=True)