#TODO contaminte xim

import os

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Sofware to contaminate a fiducial cosmology using , a dxip contaminant coming from PSF modelling')
    
    parser.add_argument('--original',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/cosmosis_pipe/2pt_sim_1110_baseline_Y3cov.fits',
                        help='File containing xip to be modified')
    parser.add_argument('--contaminant',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/abe_dxip.fits',
                        help='Path for the outputs of this code')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/alpha-beta-gamma/cosmosis_pipe/', help='location of the output of the files')
    parser.add_argument('--filename',
                        default='2pt_sim_1110_baseline_Y3cov_contaminated.fits',
                        help='Path for the outputs of this code')  
    
    args = parser.parse_args()
    return args

def main():
    import fitsio
    import itertools
    import numpy as np
    from fitsio import FITS,FITSHDR
    from astropy.io import fits
    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
        
    
    fiducialfit = args.original
    covmatrixfit_ori=fitsio.read(fiducialfit,ext=1)
    xipfit_ori=fitsio.read(fiducialfit,ext=2)

    contaminantfit = args.contaminant
    covmatrixfit_cont=fitsio.read(contaminantfit,ext=1)
    xipfit_cont=fitsio.read(contaminantfit,ext=2)
    #print(fitcontable)

    nrows = 20
    nbins=4
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin1 = (xipfit_ori['BIN1']==i)&(xipfit_ori['BIN2']==j)
        dxipbin = xipfit_cont['VALUE']
        if(len(dxipbin ) !=0 ):
            xipfit_ori['VALUE'][bin1] -=dxipbin

    hdulist = fits.open(args.original)
    #delete all xip but saving header
    oldheaders =  [hdulist[2].header]
    hdulist.pop(index=2);
 
    xiphdu = fits.BinTableHDU(xipfit_ori, name='xipt')
    hdulist.insert(2, xiphdu)
    hdulist[2].header = oldheaders[2]
    hdulist.writeto(outpath + args.filename, clobber=True)
    
if __name__ == "__main__":
    main()
