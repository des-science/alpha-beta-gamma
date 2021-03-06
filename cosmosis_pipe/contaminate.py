#TODO contaminte xim

import os

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Sofware to contaminate a fiducial cosmology using , a dxip contaminant coming from PSF modelling')
    
    parser.add_argument('--original',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/cosmosis_pipe/2pt_sim_1110_baseline_Y3cov.fits',
                        help='File containing xip to be modified')
    parser.add_argument('--contaminant',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/ab_dxi.fits',
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
    import itertools
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
    ximfit_ori=fitsio.read(fiducialfit,ext=3)

    contaminantfit = args.contaminant
    xipfit_cont=fitsio.read(contaminantfit,ext=2)
    dxipbin = xipfit_cont['VALUE']
    #print(fitcontable)

    nrows = 20
    nbins=4
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        binp = (xipfit_ori['BIN1']==i)&(xipfit_ori['BIN2']==j)
        binm = (ximfit_ori['BIN1']==i)&(ximfit_ori['BIN2']==j)
        idxbinsp =  list(itertools.compress(xrange(len(binp)),  binp))
        idxbinsm =  list(itertools.compress(xrange(len(binm)),  binm))
        if (len(idxbinsp)!=0): xipfit_ori['VALUE'][binp] -=dxipbin
        if (len(idxbinsm)!=0): ximfit_ori['VALUE'][binm] -=dxipbin

    hdulist = fits.open(args.original)
    #delete all xip but saving header
    oldheaders =  [hdulist[2].header, hdulist[3].header]
    hdulist.pop(index=2);
    hdulist.pop(index=2);
 
    xiphdu = fits.BinTableHDU(xipfit_ori)
    ximhdu = fits.BinTableHDU(ximfit_ori)
    hdulist.insert(2, xiphdu)
    hdulist.insert(3, ximhdu)
    hdulist[2].header = oldheaders[0]
    hdulist[3].header = oldheaders[1]
    hdulist.writeto(outpath + args.filename, clobber=True)
    
if __name__ == "__main__":
    main()
