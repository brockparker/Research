#! /usr/bin/env python
# use this to edit out bad regions of optical WIRO spectra that may contain atmospheric absorption 
# around 6300 and the interstellar sodium regino near 5880 Ang
# and replace with 1.0
 
import argparse
import os.path
from os import makedirs
from astropy.io import fits
import numpy as np
# need these for photometry

lam1=5885
lam2=5902
lam3=6275
lam4=6330


def main():
    parser = argparse.ArgumentParser(description='Perform photometry of images')
    parser.add_argument('filenames',metavar='filenames',nargs='+',help='List of files to photometer')
    
    args = parser.parse_args()
    
    print()
    print(args.filenames)
    print()
    
    if len(args.filenames) == 1:
        # one input filename implies that this file is a list (since why would you try to stack one file?)
        with open(args.filenames[0],'r') as f:
            filenames = [x.strip() for x in f.readlines()]
        args.filenames = filenames

    myfilenames=args.filenames
#   print(myfilenames)
#   loop over filenames
    nfiles=np.size(myfilenames)
#   print(nfiles)

    for i in range(0,int(np.size(myfilenames))) : # iterate over all images
      hdu = fits.open(myfilenames[i])
      data=hdu[0].data
      crval1=hdu[0].header['CRVAL1']
      cdelt1=hdu[0].header['CDELT1']
      lams=crval1+np.arange(0,np.size(data)) * cdelt1
      t=np.where((lams > lam1) & (lams < lam2))
      t1=np.where((lams > lam3) & (lams < lam4))
      hdu[0].data[t]=1.0
      hdu[0].data[t1]=1.0
      newfile=myfilenames[i].replace(".fits","_ED.fits")
      print(newfile)
      hdu.writeto(newfile)

if __name__ == '__main__':
    main()
