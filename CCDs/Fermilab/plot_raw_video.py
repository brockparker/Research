# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 09:37:54 2024

@author: Brock
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from astropy import units as u
from astropy.visualization import quantity_support
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from itertools import zip_longest

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def gaussian(x, a, u, s):
    return a / (np.sqrt(2 * np.pi) * s) * np.exp(-np.power((x - u) / s, 2.0) / 2)

def find_subarray(a, size, b):
    shape = a.shape[:-1] + (a.shape[-1] - size + 1, size)
    strides = a.strides + (a.strides[-1],)
    rolling_window = (np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides) == b)
    bool_indices = np.where(np.all(rolling_window, axis=1))[0]
    return bool_indices

def get_read_noise(file, amplifier):
    hdul = fits.open(file)
    overscan = 371
    gain = 1
    std_start = 750
    image = hdul[amplifier].data[1:,overscan:]/gain #371    
    header = hdul[amplifier + 1].header
    hdul.close()
        
    nsamp = int(header['NSAMP'])
        
    counts, bins = np.histogram(image.flatten(), bins=200)
    points = moving_average(bins, 2)
    
    peaks, _ = find_peaks(counts, distance=20, height=20, width=10)
    
# =============================================================================
#     for peak in peaks:
#         plt.axvline(bins[peak], color='k')
#     plt.plot(points, counts)
#     plt.show()
# =============================================================================

    if (len(peaks) != 0):
        start_peak = points[peaks[0]]
    else:
        start_peak = -2

    start_std = std_start/np.sqrt(nsamp)
        
    popt, pcov = curve_fit(gaussian, points, counts, p0=[1000000, start_peak, start_std], bounds=([0,-np.infty,0],[np.infty, np.infty, np.infty]))
    perr = np.sqrt(np.diag(pcov))

    mean = popt[1]
    standerr = popt[2]
    new_range = (mean - 5*standerr, mean + 5*standerr)
    
    ### BP Recalculate histogram centered on first peak
    fig, ax = plt.subplots(1, 1, figsize=(6,4), layout='tight')
    
    counts, bins, _ = ax.hist(image.flatten(), bins=200, range=new_range, color='slategrey')
    points = moving_average(bins, 2)
        
    popt, pcov = curve_fit(gaussian, points, counts, p0=popt, bounds=([0,-np.infty,0],[np.infty, np.infty, np.infty]))
    perr = np.sqrt(np.diag(pcov))
    
    x = np.linspace(new_range[0], new_range[1], 10000)
    ax.plot(x, gaussian(x, *popt), color='r', label = 'A={:.3f}, $\mu$={:.3f}, $\sigma$={:.3f}'.format(*popt))
    ax.set_xlabel(r'Pixel Value [e-]')
    ax.set_ylabel(r'Counts')    
    ax.set_title('Image Histogram: {}\npinit: {:.0f}, sinit: {:.0f}, psamp: {:.0f}, ssamp: {:.0f}'.format(file_name, pinit, sinit, psamp, ssamp))
    ax.legend(loc = 'upper right')
    ax.set_xlim(new_range[0], new_range[1])
    ax.tick_params(axis='both', direction='in', which='both')
    fig.tight_layout()
    plt.savefig(base_dir + '\\histogram_' + file_name + '_' + str(amplifier), dpi=250)
    plt.show()
    
    return popt[2]


plt.rc('axes', labelsize=14)
plt.rc('figure', titlesize=16)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('legend', fontsize=12)
plt.rc('xtick.major', size=6)    # size of the tick markers
plt.rc('ytick.major', size=6)    # size of the tick markers
plt.rc('xtick.minor', size=4)    # size of the tick markers
plt.rc('ytick.minor', size=4)    # size of the tick markers

amplifier = 1 # from 0 to 3

base_dir = r'C:\Users\Brock\Documents\Git\Research\CCDs\Fermilab\LTA_Data\raw_video'

pattern = '\image_lta_skip*'

hdu_files = sorted(glob.glob(base_dir + pattern + '\HDU_' + str(amplifier) + '_*.csv'))
image_files = sorted(glob.glob(base_dir + pattern + '.fz'))

quantity_support()

freq = 15*u.MHz

sinit_arr = []
pinit_arr = []
ssamp_arr = []
psamp_arr = []
std_arr = []

# TODO plot CDS sampling parmaters

for hdu_file, image_file in zip_longest(hdu_files, image_files):
    hdul = fits.open(image_file)
    image = hdul[amplifier + 1].data
    header = hdul[amplifier + 1].header
    hdul.close()
    
    
    file_name = hdu_file[hdu_file.find('image_lta_skip'):hdu_file.find(r'\HDU')]  
    
    sinit = int(header['SINIT'])
    pinit = int(header['PINIT'])
    ssamp = int(header['SSAMP'])
    psamp = int(header['PSAMP'])
    delay_integ = int(header['DELAY_INTEG_WIDTH'])
    
    read_noise = get_read_noise(image_file, amplifier)

    sinit_arr.append(sinit)
    pinit_arr.append(pinit)
    ssamp_arr.append(ssamp)
    psamp_arr.append(psamp)
    std_arr.append(read_noise)
    
    offset = 3200
    
    df = pd.read_csv(hdu_file, header=None).to_numpy()
    integ = df[offset:,2]
    waveform = df[offset:,3]/1e5

    start = [2, 0]
    ind = find_subarray(integ, 2, start) + 1
    lim = ind[0]
    length = 3000
    
    integ_crop = integ[lim:lim+length]
    waveform_crop = waveform[lim:lim+length]
    
    time = np.arange(0, length)
    time = (time/freq).to(u.us)
    

    fig, ax = plt.subplots(1, 1, figsize=(6,4), layout='tight')
    
    ax.plot(time, integ_crop, label='Integration Width', color='grey', ls=':')
    ax.plot(time, waveform_crop, label='Waveform', color='k')
    
    pinit_idx = find_subarray(integ_crop, 2, [1,0])
    
    for pinit_start_idx in pinit_idx:
        pinit_time = time[pinit_start_idx-pinit:pinit_start_idx]
        
        psamp_start_idx = pinit_start_idx - pinit
        psamp_idx = slice(psamp_start_idx - psamp,psamp_start_idx)
        
        #ax.plot(pinit_time,[1]*len(pinit_time))#,label = 'pinit')
        ax.plot(time[psamp_idx],waveform_crop[psamp_idx], label='Pedestal Sample', color='r')
        
        
    sinit_idx = find_subarray(integ_crop, 2, [0,2])

    for sinit_start_idx in sinit_idx:
        sinit_time = time[sinit_start_idx:sinit_start_idx+sinit]

        ssamp_start_idx = sinit_start_idx + sinit
        ssamp_idx = slice(ssamp_start_idx,ssamp_start_idx + ssamp)

        #ax.plot(sinit_time,[2]*len(sinit_time))#,label = 'sinit')
        ax.plot(time[ssamp_idx],waveform_crop[ssamp_idx], label='Signal Sample', color='g')

    ax.set_xlabel(r'Time [$\mu$s]')
    ax.set_ylabel(r'Pixel Value')
    ax.set_title('LTA Waveform: {}\npinit: {:.0f}, sinit: {:.0f}, psamp: {:.0f}, ssamp: {:.0f}'.format(file_name, pinit, sinit, psamp, ssamp))
    ax.tick_params(axis='both', direction='in', which='both')
    ax.set_ylim(top = 2.5)
    
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc = 'upper left')
    
    fig.tight_layout()
    # BP Invert x axis.
    
    plt.savefig(base_dir + '\\oscilloscope_' + file_name + '_' + str(amplifier), dpi=250)
    # BP Save the image.
    
    plt.show()
    
comb_arr = np.array([pinit_arr, sinit_arr, psamp_arr, ssamp_arr, std_arr]).T
    
df = pd.DataFrame(comb_arr, columns=['pinit', 'sinit', 'psamp', 'ssamp', 'std'])