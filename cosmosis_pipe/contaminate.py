import os

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Sofware to contaminate a fiducial cosmology using , a dxip contaminant coming from PSF modelling')
    
    parser.add_argument('--original',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/cosmosis_pipe/2pt_sim_1110_baseline_Y3cov.fits',
                        help='File containing xip to be modified')
    parser.add_argument('--contaminant',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/cosmosis_pipe/abg_dxip_tomo.fits',
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
        bin2 = (xipfit_cont['BIN1']==i)&(xipfit_cont['BIN2']==j)
        dxipbin = xipfit_cont['VALUE'][bin2]
        idxbins1 =  list(itertools.compress(xrange(len(bin1)),  bin1))
        idxbins2 =  list(itertools.compress(xrange(len(bin2)),  bin2))
        if (len(idxbins2)!=0):
            lowi1 =idxbins1[0]; upi1 = idxbins1[-1]
            lowi2 =idxbins2[0]; upi2 = idxbins2[-1]
            covmatrixfit_ori[lowi1:upi1,lowi1:upi1] += covmatrixfit_cont[lowi2:upi2,lowi2:upi2]
            print(np.diag(covmatrixfit_cont[lowi2:upi2,lowi2:upi2]))
        if(len(dxipbin ) !=0 ):
            xipfit_ori['VALUE'][bin1] -=dxipbin

    hdulist = fits.open(args.original)
    #delete all covariance and xip
    hdulist.pop(index=1);
    hdulist.pop(index=1);
    print(hdulist)

    covmathdu = fits.ImageHDU(covmatrixfit_ori, name='COVMAT')
    hdulist.insert(1, covmathdu)
    xiphdu = fits.BinTableHDU(xipfit_ori, name='xipt')
    hdulist.insert(2, xiphdu)
    print(hdulist)
 
    hdulist.writeto(outpath + args.filename, clobber=True)
    
if __name__ == "__main__":
    main()
