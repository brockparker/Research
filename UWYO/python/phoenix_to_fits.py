import os
os.chdir('/d/tel2/brock/Data/HST/cyc29/data/')
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os.path
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter

cdelt1=0.10 # fixed pixel grid step in Ang
sigma=7 # PHOENIX is 0.006A pixels in the optical or 0.01 in the red;
        #convolve to 0.1Ang or 16 pix FWHM = 16/2.35=7 pix  

def convert(filename):
    flux, header = fits.getdata('Orig/'+filename, header=True)
    wavs, header2 = fits.getdata('WAVE_PHOENIX-ACES-AGSS-COND-2011.fits',header=True) # wavleneght fits file
    t=np.argwhere((wavs > 3900.) & (wavs < 7600.))
    Kwavs=np.ravel(wavs[t])
    Kflux=np.ravel(flux[t])
    cKflux=gaussian_filter(Kflux,sigma) # convolve to 0.1 A per pix for optical
    crval1=Kwavs[0]                                                            
                      
                          
    lastpix=Kwavs[np.size(Kwavs)-
1]                                                                   
                         
    crpix1=1                                                                   
                      
                         
    lam=np.arange(crval1,lastpix,cdelt1) #make regular grid in
Ang                                    
                         
# add keywrods for spectra to header
    header.set('CDELT1', cdelt1)
    header.set('CRVAL1', crval1)
    header.set('CRPIX1', crpix1)
# performs spline itnterp on regular lam grid
    cs = CubicSpline(Kwavs, Kflux)
    newy=cs(lam)
# create new fits string to write
    newfile=filename.replace('.PHOENIX-ACES-AGSS-COND-2011-
HiRes.fits','_Opt.fits')
    newfile2=newfile.replace('lte','Opt/lte')
    fits.writeto(newfile2,newy,header,overwrite=True)
    print('writing:',newfile2)

def main():
    parser = argparse.ArgumentParser(description='convert PHEONIX to fits')
    parser.add_argument('filenames',metavar='filenames',nargs='+',help='List of
files to convert')
    args = parser.parse_args()
    if len(args.filenames) == 1:
        # one input filename implies that this file is a list (since why would
you try to stack one fil
e?)
        with open(args.filenames[0],'r') as f:
            filenames = [x.strip() for x in f.readlines()]
        args.filenames = filenames

    myfilenames=args.filenames
 # loop over filenames
    nfiles=np.size(myfilenames)
    print('Read ',nfiles,' FITS files from ',args.filenames)

    for i in range(0,int(np.size(myfilenames))) : # iterate over all images
       convert(myfilenames[i])
          
if __name__ == '__main__':
    main()