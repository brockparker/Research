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

folder = '/home/GitHub/Research/CCDs/Sophia_QE/'

dates = ['','','','20240724']

wl = 300

prefix = 'warm'

science_file = r'D:\Sophia\20240724\QE\test_science_{:.0f}.fits'.format(wl)
dark_file = r'D:\Sophia\20240724\QE\test_dark_{:.0f}.fits'.format(wl)
bias_file = r'D:\Sophia\20240724\QE\test_bias_{:.0f}.fits'.format(wl)

science = fits.getdata(science_file)[0].astype('int32')
dark = fits.getdata(dark_file)[0].astype('int32')
bias = fits.getdata(bias_file)[0].astype('int32')

reduced = (science - bias) - (dark - bias)

plt.imshow(bias, vmin = 1000, vmax = 25000, cmap='inferno')
plt.title('Bias')
plt.colorbar()
plt.savefig(r'D:\Sophia\20240724\QE\bias_{}_{:.0f}.png'.format(prefix, wl), dpi = 250)
plt.show()

plt.imshow(science, vmin = 10000, vmax = 60000, cmap='inferno')
plt.title('Raw Science')
plt.colorbar()
plt.savefig(r'D:\Sophia\20240724\QE\science_raw_{}_{:.0f}.png'.format(prefix, wl), dpi = 250)
plt.show()

plt.imshow(dark, vmin = 10000, vmax = 30000, cmap='inferno')
plt.title('Raw Dark')
plt.colorbar()
plt.savefig(r'D:\Sophia\20240724\QE\dark_raw_{}_{:.0f}.png'.format(prefix, wl), dpi = 250)
plt.show()

plt.imshow(dark - bias, vmin = 6000, vmax = 8000, cmap='inferno')
plt.title('Bias Subtracted Dark')
plt.colorbar()
plt.savefig(r'D:\Sophia\20240724\QE\dark_sub_bias_{}_{:.0f}.png'.format(prefix, wl), dpi = 250)
plt.show()


plt.imshow(reduced, norm=LogNorm(vmin = 1000, vmax = 30000), cmap='inferno')
plt.title('Reduced Science')
plt.colorbar()
plt.savefig(r'D:\Sophia\20240724\QE\science_reduced_{}_{:.0f}.png'.format(prefix, wl), dpi = 250)
plt.show()