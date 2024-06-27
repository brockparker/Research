'''
Skipper Read Noise Modelling 
by Brock Parker

Updated 5/20/2024
'''

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c
from astropy.visualization import quantity_support
import pandas as pd
from scipy.integrate import simpson
from scipy.optimize import curve_fit
from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
from astropy.modeling import Fittable1DModel, Parameter
from scipy.signal import savgol_filter
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
# BP Import needed packages.

machine = 'Laptop'
# BP Define what machine code is being run on to save files. For internal use only.

if machine == 'Linux':
    path = '/home/baparker/GitHub'
elif machine == 'Laptop':
    path = 'C:/Users/Brock/Documents/Git'
elif machine == 'Desktop':
    path = 'F:/Github'
# BP Define file save path based on machine.
    
plt.rc('axes', labelsize=14)
plt.rc('figure', titlesize=30)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('legend', fontsize=12)
# BP Plot stylization parameters.

plot_colors = ['#293462', '#1CD6CE', '#FEDB39', '#D61C4E', '#0b132b']

@u.quantity_input(f=u.Hz, tau=u.s, delay=u.s)
def transfer(f, tau, n, delay = 0*u.s):
    '''
    Function defining the transfer function for correlated double sampling.
    Valid for both single samples and multiple samples at varying pedestal sampling times.
    
    f : astropy.units.Quantity
        Array of input frequencies as an astropy Quantity. This is the value the function is evaluated at, and represents the specific frequency contribution.
    
    tau : astropy.units.Quantity
        The pedestal integration time, i.e. how long it takes to move the charge on and off of the pedestal. The individual pixel time is thus half this.
        
    n : int
        The number of times a sample is taken. For Skipper CCDs, this is analagous to the number of skips.
    '''
    
    tau = tau / 2

    part1 = 4 * (np.sin((np.pi * f * tau).to(u.rad, equivalencies=u.dimensionless_angles())))**2 / (2 * np.pi * f * tau * n).to(u.rad, equivalencies=u.dimensionless_angles())
    part2 = abs( np.sin((2*np.pi * f * tau * n).to(u.rad, equivalencies=u.dimensionless_angles())) / np.sin((2*np.pi * f * tau).to(u.rad, equivalencies=u.dimensionless_angles())))

    transfer_value = part1*part2

    return transfer_value
# BP Define function to return the transfer function. Describes how much weight is given to noise at each frequency.

@u.quantity_input(ewn=u.V/u.Hz**0.5, f1=u.Hz, f2=u.Hz, tau=u.s, delay=u.s, conversion=u.V/u.electron)
def total_noise_model(n, ewn, tau, f1, f2, conversion, delay=0*u.s):
    '''
    Theoretical total noise model as a function of noise spectrum parameters. Valid for noise model as defined by the 1/f and 1/f^2 noise corner and white noise level.
    Defined in Greffe & Smith 2024 Poster.
    Valid for both single and mutiple samples.
    
    n : int
        The number of times a sample is taken. For Skipper CCDs, this is analagous to the number of skips.
    
    ewn : astropy.units.Quantity
        White noise level in V/Hz^1/2. Should be the value the noise spectrum levels off at.
    
    tau : astropy.units.Quantity
        The pedestal integration time, i.e. how long it takes to move the charge on and off of the pedestal. The individual pixel time is thus half this.
    
    f1 : astropy.units.Quantity
        1/f noise corner frequency. Where the noise spectrum turns over.
        
    f2 : astropy.units.Quantity
        1/f^2 noise corner frequency. Where the noise spectrum turns over.
        
    conversion : astropy.units.Quantity
        Conversion factor to convert between number of electrons and volts.
    '''
    tau = tau / 2
    
    total_noise = ewn * np.sqrt(1/tau + 4*f1*np.log(2) + 4/3 * f2**2 * np.pi**2 * tau) / np.sqrt(n)

    return total_noise / conversion

@u.quantity_input(f=u.Hz, tau=u.s, noise_density=u.V/u.Hz**0.5, conversion=u.V/u.electron)
def total_noise_integral(f, tau, noise_density, conversion, n):
    '''
    The total noise calculated by convolving the noise density spectrum with the transfer function.
    All values must be in MKS.
    
    f : astropy.units.Quantity
        Array of input frequencies as an astropy Quantity. This is the value the function is evaluated at, and represents the specific frequency contribution.
    
    tau : astropy.units.Quantity
        The pedestal integration time, i.e. how long it takes to move the charge on and off of the pedestal. The individual pixel time is thus half this.
        
    noise_density : astropy.units.Quantity
        The noise density spectrum from the on-chip MOSFET. In units of V/Hz^1/2     
        
    conversion : astropy.units.Quantity
        Conversion factor to convert between number of electrons and volts.
    '''

    try:
        variance = []
        for i in range(len(tau)):
            value = simpson((noise_density * transfer(f, tau[i], n))**2, f)
            variance.append(value)

    except:
        variance = simpson((noise_density * transfer(f, tau, n))**2, f)
        
    sigma = np.sqrt(variance)
    return sigma * u.V / conversion

@u.quantity_input(f=u.Hz, f1=u.Hz, f2=u.Hz, conversion=u.V/u.electron, ewn=u.V/u.Hz**0.5)
def noise_spectrum_model(f, ewn, f1, f2, conversion):
    '''
    Theoretical MOSFET noise spectrum function to allow extrapolation from measured spectrum by best fit.
    
    f : astropy.units.Quantity
        Array of input frequencies as an astropy Quantity. This is the value the function is evaluated at, and represents the specific frequency contribution.
    
    ewn : astropy.units.Quantity
        White noise level in V/Hz^1/2. Should be the value the noise spectrum levels off at.
    
    f1 : astropy.units.Quantity
        1/f noise corner frequency. Where the noise spectrum turns over.
        
    f2 : astropy.units.Quantity
        1/f^2 noise corner frequency. Where the noise spectrum turns over.
        
    conversion : astropy.units.Quantity
        Conversion factor to convert between number of electrons and volts.
    '''

    return ewn * (1/(1 + (f/f2)**2) + 1/(f/f1)) / conversion

quantity_support()
# BP Allow plotting with units.

fig, ax = plt.subplots(1, 1, figsize=(5,4), layout='tight')
# BP Initialize plots.

frequencies = np.logspace(0,10,100000) * u.Hz
taus = np.logspace(-10, 0, 100000) * u.s 
# BP Create array of frequencies and pedestal times.

tau_10 = 1/(10*u.kHz)
tau_100 = 1/(100*u.kHz)
tau_1000 = 1/(1*u.MHz)
# BP Define three discrete pedestal times.

ax.plot(frequencies, transfer(frequencies, tau_10, 1), label = '{:.0f}'.format((1/tau_10).to(u.kHz)), color = plot_colors[0], ls='-')
ax.plot(frequencies, transfer(frequencies, tau_100, 1), label = '{:.0f}'.format((1/tau_100).to(u.kHz)), color = plot_colors[1], ls='-')
ax.plot(frequencies, transfer(frequencies, tau_1000, 1), label = '{:.0f}'.format((1/tau_1000).to(u.kHz)), color = plot_colors[2], ls='-')
# BP Plot transfer functions.

ax.set_xlabel(r'Frequency [Hz]')
ax.set_ylabel(r'Transfer Function')
ax.legend(loc = 'best')
ax.set_xlim(1e3, 0.2e7)
ax.set_xscale('log')
ax.set_ylim(0, 1.75)
ax.tick_params(axis='both', direction='in', which='both')
fig.tight_layout()
# BP Sylization parameters.

plt.savefig(path + '/Research/CCDs/Skipper_Noise/cds_transfer_func.png', dpi=250)
plt.show()
# BP Save figure.

# =============================================================================
# stop
# =============================================================================

fig, ax = plt.subplots(1, 1, figsize=(5,4), layout='tight')
# BP Initialize plots.

tau_200 = 1/(2e5*u.Hz)
n_skip_1 = 2
n_skip_2 = 32
# BP Define pedestal time and number of skips.

ax.plot(frequencies, transfer(frequencies, tau_200, 1), color = plot_colors[0], ls='-', label='1 Skip')
ax.plot(frequencies, transfer(frequencies, tau_200, n_skip_1), color = plot_colors[1], ls='-', label='2 Skips')
ax.plot(frequencies, transfer(frequencies, tau_200, n_skip_2), color = plot_colors[2], ls='-', label='32 Skips')
# BP Plot transfer function of different number of skips.

ax.set_xlabel(r'Frequency [Hz]')
ax.set_ylabel(r'Transfer Function')
ax.legend(loc = 'best')
ax.set_xlim(1e3, 1e7)
ax.set_xscale('log')
ax.set_ylim(0, 1.75)
ax.tick_params(axis='both', direction='in', which='both')
fig.tight_layout()
# BP Stylization parameters.

plt.savefig(path + '/Research/CCDs/Skipper_Noise/skipper_cds_transfer_func.png', dpi=250)
plt.show()
# BP Save figure.

# =============================================================================
# stop
# =============================================================================

fig, ax = plt.subplots(1, 1, figsize=(5,4), layout='tight')
ax2 = ax.twiny()
# BP Initialize plots.

wn_floor = 2e-8 * u.V / ((u.Hz)**(1/2))
# BP Define white noise floor.

nc_1f = 100 * u.kHz
nc_1f2 = 100 * u.kHz#2e6 * u.Hz # Nominally 100 * u.kHz
# BP Define noise corners.

n_skips = 1024
# BP Define number of skips.

n_pix_x = 4000
n_pix_y = 2000
n_pix = n_pix_x * n_pix_y
# BP Define number of pixels.

n_amp = 128
# BP Define number of amplifiers.

conversion = 2.5e-6 * u.V/u.electron
# BP Steve recommended 5e-6
# BP LTA measured at 2.5e6

ax.plot(taus, total_noise_model(1, wn_floor, taus, nc_1f, nc_1f2, conversion), label='Traditional CCD', color=plot_colors[0])
ax.plot(taus * n_skips, total_noise_model(n_skips, wn_floor, taus, nc_1f, nc_1f2, conversion), label='{:.0f} Skips'.format(n_skips), color=plot_colors[0], ls='--')
ax.axhline(0.15, color=plot_colors[3], ls=':', label='0.15 e$^-$ threshold')
# BP Plot total integrated noise using theoretical model.

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'Pixel Time [s]')
ax.set_ylabel(r'Integrated Noise [e-]')
ax.legend(loc = 'best')
ax.set_xlim(1e-6, 1e-1)
ax.set_ylim(5e-3, 1e2)
ax.tick_params(axis='both', direction='in', which='both')
# BP Pedestal time axes sytlization.

ax2.set_xscale('log')
ax2.set_xlabel(r'Readout Time [s]')
ax2.set_xlim(ax.get_xlim()[0] * n_pix / n_amp, ax.get_xlim()[1] * n_pix / n_amp)
ax2.tick_params(axis='both', direction='in', which='both')
fig.tight_layout()
# Calculate corresponding read time for secondary axes. Stylization parameters.

plt.savefig(path + '/Research/CCDs/Skipper_Noise/model_integrated_noise_128amp_4kx2k.png', dpi=250)
plt.show()
# BP Save figure.

# =============================================================================
# stop
# =============================================================================

noise_spectrum = pd.read_csv(path + '/Research/CCDs/Skipper_Noise/STANoiseSpectrum.csv', header=None)
lta_noise_spectrum = pd.read_csv(path + '/Research/CCDs/Skipper_Noise/lta_noise_spectrum_3.csv', header=None)
# BP Read in read-out amplifier noise spectrums.

noise_freqs = np.array(noise_spectrum[0]) * u.Hz 
noise_density = np.array(noise_spectrum[1]) * u.V / (u.Hz**0.5)
# BP Extract noise frequencies and densities from data.

lta_noise_freqs = np.array(lta_noise_spectrum[0]) * u.Hz * 200
lta_noise_density = np.array(lta_noise_spectrum[1]) * 1 # Preamplifier gain of 10
lta_noise_density = (savgol_filter(lta_noise_density, 50, 3) * u.nV / (u.Hz**0.5)).to(u.V / (u.Hz**0.5))
# BP Extract noise frequencies and densities from data.

opt_speed = 1/(1000*u.kHz)
cds_transfer = transfer(frequencies, opt_speed, 1)
cds_transfer_32 = transfer(frequencies, opt_speed, n_skip_2)
# BP Define optimum readout speed and create CDS transfer function at that frequency.

fig, ax = plt.subplots(1, 1, figsize=(6,4), layout='tight')
ax2 = ax.twinx()
# BP Create plot axes.

ax2.plot(frequencies, cds_transfer, color = plot_colors[0], ls='--', label='1 Skip')
ax2.plot(frequencies, cds_transfer_32, color = plot_colors[3], ls='--', label='32 Skips')
ax.plot(noise_freqs, noise_density, color = plot_colors[2], ls='-', label='STA')
ax.plot(lta_noise_freqs, lta_noise_density, color = plot_colors[1], ls='-', label='LTA')
ax.plot(0, 0, color = plot_colors[0], ls='--', label='1 Skip')
ax.plot(0, 0, color = plot_colors[3], ls='--', label='32 Skips')
# BP plot noise spectrum and CDS function.

ax.set_xlabel(r'Frequency [Hz]')
ax2.set_ylabel(r'Transfer Function')
ax.set_ylabel(r'Noise Density [V Hz$^{-1/2}$]')
ax.legend(loc = 'lower left')
ax.set_xlim(2.1e3, 1e7)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1e-10, 1e-5)
ax2.set_ylim(0, 2)
ax.tick_params(axis='both', direction='in', which='both')
fig.tight_layout()
# BP Stylization parameters.

plt.savefig(path + '/Research/CCDs/Skipper_Noise/measured_noise_spectrums.png', dpi=250)
plt.show()
# BP Save figure.

# =============================================================================
# stop
# =============================================================================

fig, ax = plt.subplots(1, 1, figsize=(6,4.5), layout='tight')
ax2 = ax.twiny()
# BP Create plotting axes.

n_skips_64 = 64

array_len = 20
# BP Number of points in plotted noise_spectrum.

taus = 1 / (np.logspace(0, 7, array_len) * u.Hz)
plot_taus = 1 / (np.logspace(3.5, 6.5, array_len) * u.Hz)
# BP Create plotting variables - different lenght than calculation range.

output_noise_sta = total_noise_integral(noise_freqs, plot_taus, noise_density, conversion, n=1)
output_noise_lta = total_noise_integral(lta_noise_freqs, plot_taus, lta_noise_density, conversion, n=1)
# BP Calculate total integrated noise from measured noise_spectrums.

ax.plot(taus, total_noise_model(1, wn_floor, taus, nc_1f, nc_1f2, conversion), color=plot_colors[0], label='Analytical Model')
ax.plot(plot_taus, output_noise_sta, color=plot_colors[2], ls='-', label='STA', marker='x')
ax.plot(plot_taus, output_noise_lta, color=plot_colors[1], ls='-', label='LTA', marker='x')
ax.plot(plot_taus * n_skips_64, output_noise_sta / np.sqrt(n_skips_64), color=plot_colors[2], ls='--', marker='')
ax.plot(plot_taus * n_skips_64, output_noise_lta / np.sqrt(n_skips_64), color=plot_colors[1], ls='--', marker='')
ax.plot(plot_taus * n_skips, output_noise_sta / np.sqrt(n_skips), color=plot_colors[2], ls=':', marker='')
ax.plot(plot_taus * n_skips, output_noise_lta / np.sqrt(n_skips), color=plot_colors[1], ls=':', marker='')
ax.axhline(0.15, color=plot_colors[3], ls=':', label='0.15 e$^-$ threshold')
ax.axhline(0.001, color=plot_colors[0], ls='--', label='64 Skips')
ax.axhline(0.001, color=plot_colors[0], ls=':', label='1024 Skips')
# BP Plot measured integrated total noise compared to

# =============================================================================
# ax2.plot(1/taus, [0]*array_len, color='green', ls='-')
# =============================================================================
# ax.axvline(1/(200 * u.kHz), label='200 kHz Optimum', color='k', ls=':')
# ax.axvline(1/(100 * u.kHz), label='100 kHz Optimum', color='k', ls=':')
# ax.axhline(0.15, color='red', ls='--', label='0.15 e$^-$ threshold')

ax.set_xlabel(r'Pixel Time [s]')
ax.set_ylabel(r'Integrated Noise [e-]')
ax.legend(loc = 'best')
ax.set_xlim(1e-6, 1e-1)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(8e-2, 3e1)
ax.tick_params(axis='both', direction='in', which='both')

ax2.set_xscale('log')
ax2.set_xlabel(r'Readout Time [s]')
ax2.set_xlim(ax.get_xlim()[0] * n_pix / n_amp, ax.get_xlim()[1] * n_pix / n_amp)
ax2.tick_params(axis='both', direction='in', which='both')
fig.tight_layout()
# BP Stylization parameters.

plt.savefig(path + '/Research/CCDs/Skipper_Noise/measured_integrated_noise.png', dpi=250)
plt.show()
# Bp Save figure.

stop

wn_floor = 2e-8 * u.V / ((u.Hz)**(1/2))

nc_1f = 100 * u.kHz
nc_1f2 = 100 * u.kHz#2e6 * u.Hz # Nominally 100 * u.kHz

best_guess = [1.9e-8, 1e5, 0.95e7]

popt, pcov = curve_fit(noise_spectrum_model, noise_freqs, noise_density, p0=best_guess, bounds=([0, 0, 0],[np.infty,np.infty,np.infty]), sigma = np.sqrt(noise_density), absolute_sigma=False)
cerr = np.sqrt(np.diag(pcov))

sig_below = noise_spectrum_model(frequencies, popt[0]-cerr[0],popt[1]-cerr[1],popt[2]-cerr[2], conversion)
sig_above = noise_spectrum_model(frequencies, popt[0]+cerr[0],popt[1]+cerr[1],popt[2]+cerr[2], conversion)

# =============================================================================
# z = np.polyfit(noise_freqs.value, noise_density.value, 7)
# p = np.poly1d(z)
# =============================================================================




# =============================================================================
# model = noise_model(*best_guess)  
# fitter = fitting.LevMarLSQFitter()
# noise_fit = fitter(model, noise_freqs, noise_density)
# =============================================================================



fig, ax = plt.subplots(1, 1, figsize=(6,4.5), layout='tight')

popt2 = [1.9e-8, 1e5, 0.9e7]

ax.plot(frequencies, noise_spectrum_model(frequencies, *popt, conversion), label='Total Noise \ne$_{{wn}}$={:.2e}$\pm${:.2e}\nf$_{{1nc}}$={:.2e}$\pm${:.2e}\n$f_{{2nc}}=${:.2e}$\pm${:.2e}'.format(popt[0],cerr[0],popt[1],cerr[1],popt[2],cerr[2]), color='purple')
ax.plot(frequencies, noise_spectrum_model(frequencies, *popt2, conversion), label='Theoretical', color='green')
#ax.plot(freqs.value, p(freqs.value), color='orange')
ax.plot(noise_freqs, noise_density, color='red', ls='-')

ax.fill_between(frequencies, sig_below, sig_above, facecolor='purple', alpha=0.25)

ax.set_xlabel(r'Frequency [Hz]')
ax.set_ylabel(r'Noise Density [nV Hz$^{-1/2}$]')
ax.legend(loc = 'upper right')
ax.set_xlim(1e2, 1e8)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1e-9, 1e-5)
ax.axhline(wn_floor)

ax.tick_params(axis='both', direction='in', which='both')

fig.tight_layout()

plt.savefig(path + '/Research/CCDs/Skipper_Noise/best_fit_noise_spectrum.png', dpi=250)
plt.show()






array_len = 500

model_noise_density = noise_spectrum_model(frequencies, *popt2)

output_noise = total_noise_integral(taus, noise_density, noise_freqs, conversion)
model_output_noise = total_noise_integral(frequencies, taus, model_noise_density, conversion)

fig, ax = plt.subplots(1, 1, figsize=(6,4.5), layout='tight')
# ax3 = ax.twiny()

ax.plot(taus, output_noise / conversion, color='blue', ls='', label='Real Noise', marker='*')
ax.plot(taus, model_output_noise / conversion, color='red', ls='--', label='Model Noise')
ax.plot(taus, total_noise_model(wn_floor, taus, nc_1f, nc_1f2) / conversion, color='black', label='Analytical Model')

# ax3.plot(taus, [0]*array_len, color='green', ls='-')

# ax.axvline(1/(200 * u.kHz), label='200 kHz Optimum', color='k', ls=':')
# ax.axhline(0.15, color='red', ls='--', label='0.15 e$^-$ threshold')

ax.set_xlabel(r'Pixel Time [s]')
ax.set_ylabel(r'Integrated Noise [e-]')
# ax.set_xlabel(r'Frequency [Hz]')
ax.legend(loc = 'best')
ax.set_xlim(1e-10, 1e-2)
ax.set_xscale('log')
ax.set_yscale('log')
# ax3.set_xscale('log')
ax.set_ylim(1e0, 1e2)

ax.tick_params(axis='both', direction='in', which='both')
 
fig.tight_layout()

plt.savefig(path + '/Research/CCDs/ccd_integrated_noise_data.png', dpi=250)
plt.show()


## Fit by CDS function



