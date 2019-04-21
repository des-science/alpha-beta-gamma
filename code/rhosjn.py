import numpy as np
import os
import kmeans_radec
import matplotlib.pyplot as plt

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    
    parser.add_argument('--piff_cat',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3a1-v29',
                        help='Full Path to the Only stars Piff catalog')
    parser.add_argument('--exps_file',
                        default='/home/dfa/sobreira/alsina/DESWL/psf/ally3.grizY',
                        #default='/home/dfa/sobreira/alsina/DESWL/psf/testexp', 
                        help='list of exposures (in lieu of separate exps)')
    parser.add_argument('--bands', default='riz', type=str,
                         help='Limit to the given bands')
    parser.add_argument('--bandcombo', default=False,
                        action='store_const', const=True,
                        help='run rho2 for all combination of bands, if false run particular combination defined in band')
    parser.add_argument('--use_reserved', default=True,
                        action='store_const', const=True,
                        help='just use the objects with the RESERVED flag')
    parser.add_argument('--frac', default=1., type=float,
                        help='Choose a random fraction of the input stars')
    parser.add_argument('--mod', default=True,
                        action='store_const', const=True,
                        help='If true it substracts the mean to each field before calculate correlations')
    parser.add_argument('--sn', default=True,
                        action='store_const', const=True,
                        help='If true multiply by 2 the variances of all correlations. Shape-noise error.')
    parser.add_argument('--xim', default=False,
                        action='store_const', const=True,
                        help='trecorr return xim instead of xip')
    parser.add_argument('--obs', default=False,
                        action='store_const', const=True,
                        help='Use e_obs instead of e_piff to calculate modified rho stats')
    parser.add_argument('--outpath', default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/',
                        help='location of the output of the files')
    parser.add_argument('--filename', default='allrhosp_jk.fits', help='Name of the fit file where info of dxip will be saved ')
    
    
    args = parser.parse_args()

    return args

def jk_kmeans(ra,dec,njk,plot=False):
	'''
	Function that takes RA and Dec from a given catalog and computes JK regions using kmeans_radec module by Erin Sheldon.

	Parameters
	----------
        ra, dec : numpy arrays of RA and Dec. len(ra) = len(dec) = number of galaxies.
	njk : number of JK regions.

	Returns
	-------
        jk = JK region for each galaxy: integer ranging from 0 to njk-1. It is numpy array with the same length as ra and dec. 

	'''
        radec = np.zeros((len(ra),2))
        radec[:,0] = ra
        radec[:,1] = dec
        km = kmeans_radec.kmeans_sample(radec,njk,maxiter=500,tol=1e-05)
        jk = km.find_nearest(radec)
        if not km.converged:
                print 'k means did not converge'
        if plot:
                plt.figure()
                plt.scatter(ra,dec,c=jk,lw=0,cmap='Paired',rasterized=True)
                plt.xlabel(r'RA',fontsize=12)
                plt.ylabel(r'Dec',fontsize=12)
                plt.savefig('jk_kmeans.pdf')
        return jk

def write_fit(data, names, filename):
    import fitsio
    from fitsio import FITS,FITSHDR
    import os.path
    if not os.path.isfile(filename):
        fits = FITS(filename,'rw')
        fits.write(data, names=names, clobber=False)
        print("Writing file: ", filename)
    else:
        fits = FITS(filename,'rw')
        fits[-1].append(data)
        print("Apending File: ",  filename)
def main():
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    #sys.path.insert(0, '/global/cscratch1/sd/alsina/alpha-beta-gamma/code/src')

    import numpy as np
    from read_psf_cats import read_data,  toList
    from run_rho import measure_rho
    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
        

    #STATISTIC USING ONLY RESERVED STARS
    keys = ['ra', 'dec','obs_e1', 'obs_e2', 'obs_T',
            'piff_e1', 'piff_e2', 'piff_T',  'mag']
 
    exps = toList(args.exps_file)
    data, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved)
    print("Objects",  len(data))
    data = data[data['mag']<20]
    print("Objects with magnitude <20",  len(data))
    
    names=['JKR', 'ANGBIN','THETA', 'RHO0', 'VAR_RHO0', 'RHO1', 'VAR_RHO1', 'RHO2', 'VAR_RHO2', 'RHO3', 'VAR_RHO3', 'RHO4', 'VAR_RHO4', 'RHO5', 'VAR_RHO5']
    forms = ['i4', 'i4', 'f8',  'f8',  'f8', 'f8',  'f8',  'f8', 'f8',  'f8',  'f8', 'f8',  'f8',  'f8', 'f8']
    dtype = dict(names = names, formats=forms)
    nrows = 20
    outdata = np.recarray((nrows, ), dtype=dtype)

    njk = 100
    jkindexes = jk_kmeans(data['ra'], data['dec'] ,njk,plot=False)
    #print (jkindexes)
    
    for jkidx in range(njk):
        print("running jackkniffe region",  jkidx)
        rho0, rho1, rho2, rho3, rho4, rho5 = measure_rho(data[ jkindexes!=jkidx ],  mod=args.mod, obs=args.obs)
        jkrarr =  np.array([jkidx]*nrows)
        angarr = np.arange(nrows)
        thetaarr = np.exp(rho0.meanlogr)
        if args.xim:
            rho0arr = rho0.xim; rho1arr = rho1.xim; rho2arr = rho2.xim; rho3arr = rho3.xim; rho4arr = rho4.xim; rho5arr = rho5.xim
        else:
            rho0arr = rho0.xip; rho1arr = rho1.xip; rho2arr = rho2.xip; rho3arr = rho3.xip; rho4arr = rho4.xip; rho5arr = rho5.xip
        varrho0arr = 2*rho0.varxi; varrho1arr = 2*rho1.varxi; varrho2arr = 2*rho2.varxi;
        varrho3arr = 2*rho3.varxi; varrho4arr = 2*rho4.varxi; varrho5arr = 2*rho5.varxi;
        array_list = [jkrarr, angarr, thetaarr, rho0arr, varrho0arr,
                      rho1arr, varrho1arr,rho2arr, varrho2arr,rho3arr,
                      varrho3arr,rho4arr, varrho4arr,rho5arr,
                      varrho5arr]
        for array, name in zip(array_list, names): outdata[name] = array
        write_fit(outdata, names, outpath + args.filename)
    
if __name__ == "__main__":
    main()
