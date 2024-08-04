# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 09:55:01 2024

@author: Brock
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import glob
from astropy.io import fits
from scipy.signal import find_peaks

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def gaussian(x, a, u, s):
    return a / (np.sqrt(2 * np.pi) * s) * np.exp(-np.power((x - u) / s, 2.0) / 2)

plt.rc('axes', labelsize=14)
plt.rc('figure', titlesize=30)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('legend', fontsize=12)
plt.rc('xtick.major', size=6)    # size of the tick markers
plt.rc('ytick.major', size=6)    # size of the tick markers
plt.rc('xtick.minor', size=4)    # size of the tick markers
plt.rc('ytick.minor', size=4)    # size of the tick markers

path = r'C:\Users\Brock\Documents\Git\Research\CCDs\Fermilab\LTA_Data\read_noise'

#files = sorted(glob.glob(path + '\proc_image_lta_reverse_??.fits'))
files = sorted(glob.glob(path + '\proc_image_lta_??.fits'))

amplifier = 1 # only amplifier 2 is working, from 1 to 4
overscan = 371 # pixel where overscan starts, TODO define better
gain = 150 # approximate gain

nsamp = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
#nsamp = [1, 2, 3, 4, 8, 12, 16, 24, 32, 45, 64, 100, 128, 200, 256, 512, 1024, 2048, 4096]
std, std_err = [], []

std_start = 10

for i, file in enumerate(files):
    file_name = file[file.find('read_noise')+11:-5]
    
    hdul = fits.open(file)
    #image = hdul[amplifier].data[1:,1:]/gain #ignore first column and row, known transient
    image = hdul[amplifier].data[1:,overscan:]/gain #371
    hdul.close()
    
    image[image > 100] = np.nan

    fig, ax = plt.subplots(1, 1, figsize=(6,4), layout='tight')

    counts, bins, _ = ax.hist(image.flatten(), bins=250, color='slategrey', label='Physical Histogram')
    points = moving_average(bins, 2)
    
    peaks, _ = find_peaks(counts, distance=10, height=10, width=3)

    #for peak in peaks:
        #plt.axvline(bins[peak], color='k')
        
    if (len(peaks) != 0):
        start_peak = points[peaks[0]]
    else:
        start_peak = -2
    start_std = std_start/np.sqrt(nsamp[i])
        
    popt, pcov = curve_fit(gaussian, points, counts, p0=[100, start_peak, start_std])#, bounds=([0,-np.infty,0],[np.infty, np.infty, np.infty]))
    perr = np.sqrt(np.diag(pcov))
    
    print(popt)
    #if file_name == 'proc_image_lta_17':
    #    stop
    x = np.linspace(image.min(), image.max(), 100000)
    ax.plot(x, gaussian(x, *popt), color='r', label = 'Fit Gaussian')
    ax.set_xlabel(r'Pixel Value [e-]')
    ax.set_ylabel(r'Counts')
    ax.set_title('Image Histogram')
    ax.legend(loc = 'upper right')
    #ax.set_xlim(-1000, 2000)
    ax.tick_params(axis='both', direction='in', which='both')
    fig.tight_layout()
    #plt.savefig(path + '\' + '_' + str(amplifier), dpi=250)
    plt.show()
    
    mean = popt[1]
    standerr = popt[2]
    new_range = (mean - 5*standerr, mean + 5*standerr)
    
    ### BP Recalculate histogram centered on first peak
    fig, ax = plt.subplots(1, 1, figsize=(6,4), layout='tight')

    counts, bins, _ = ax.hist(image.flatten(), bins=150, range=new_range, color='slategrey', label='Physical Histogram')
    points = moving_average(bins, 2)
        
    popt, pcov = curve_fit(gaussian, points, counts, p0=popt)#, bounds=([0,-np.infty,0],[np.infty, np.infty, np.infty]))
    perr = np.sqrt(np.diag(pcov))
    
    x = np.linspace(new_range[0], new_range[1], 10000)
    ax.plot(x, gaussian(x, *popt), color='r', label = 'Fit Gaussian')
    ax.set_xlabel(r'Pixel Value [e-]')
    ax.set_ylabel(r'Counts')
    ax.set_title('Image Histogram')
    ax.legend(loc = 'upper right')
    ax.set_xlim(new_range[0], new_range[1])
    ax.tick_params(axis='both', direction='in', which='both')
    fig.tight_layout()
    #plt.savefig(path + '\\')
    plt.show()
    
    std.append(popt[2])
    std_err.append(perr[2])
    

        
std = np.array(std)
std_err = np.array(std_err) * 3 # 3 sigma error bars
    
fig, ax = plt.subplots(1, 1, figsize=(6,4), layout='tight')

idx = 0 # start at NSAMP = 2 since the first image was not processed

single_sample_std = std[idx] * np.sqrt(nsamp[idx])

ax.plot(nsamp, single_sample_std/np.sqrt(nsamp), color='k', label='Theoretical', ls='--')
ax.plot(nsamp, std, color='r', label='Measured')

ax.fill_between(nsamp, std - std_err, std + std_err, color='r', alpha=0.25, edgecolor=None)

ax.set_xlabel(r'Number of Samples')
ax.set_ylabel(r'Noise [e-]')
ax.set_title('Noise vs Skips Overscan')
ax.legend(loc = 'upper right')
ax.set_xscale('log')
ax.set_yscale('log')
ax.tick_params(axis='both', direction='in', which='both')
fig.tight_layout()
#plt.savefig(path + '\\noise_curve_reversed_' + str(amplifier) + '.png', dpi=250)
plt.show()