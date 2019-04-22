import numpy as np
import os
import kmeans_radec
import matplotlib.pyplot as plt

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    
    parser.add_argument('--metacal_cat',
                        #default='/home2/dfa/sobreira/alsina/catalogs/y3_master/Y3_mastercat_v2_6_20_18_subsampled.h5',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3_master/Y3fullmaster/Y3_mastercat_v2_6_20_18.h5', 
                        help='Full Path to the Metacalibration catalog')
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
    parser.add_argument('--outpath', default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/',
                        help='location of the output of the files')
    parser.add_argument('--tomo', default=False,
                        action='store_const', const=True,
                        help='Run all tomographic correlations')
    parser.add_argument('--nz_source',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3_master/nz_source_zbin.h5',
                        help='Full Path to the Only stars Piff catalog')
    parser.add_argument('--filename', default='alltausp_jk.fits', help='Name of the fit file where info of dxip will be saved ')
    
    
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
        km = kmeans_radec.kmeans_sample(radec,njk,maxiter=100,tol=1e-01)
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
    from read_psf_cats import read_data, toList, read_h5
    from run_rho import  measure_tau
    import h5py as h
    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
        

  
    #Reading Mike stars catalog
    keys = ['ra', 'dec','obs_e1', 'obs_e2', 'obs_T',
            'piff_e1', 'piff_e2', 'piff_T', 'mag']
 
    exps = toList(args.exps_file)
    data_stars, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved)
    print("Objects",  len(data_stars))
    data_stars = data_stars[data_stars['mag']<20]
    print("Objects with magnitude <20",  len(data_stars))
    
    
    if(args.tomo):
        #Make directory where the ouput data will be
        ipath =  os.path.join(args.outpath, 'tomo_taus' )
        outpath = os.path.expanduser(ipath)
        try:
            if not os.path.exists(outpath):
                os.makedirs(outpath)
        except OSError:
            if not os.path.exists(outpath): raise
        print('Starting Tomography!')
        galkeys = ['ra','dec','e_1','e_2','R11','R22']
        data_gal =  read_h5(args.metacal_cat, 'catalog/metacal/unsheared',  galkeys )
        print("Total objects in catalog:", len(data_gal))
        dgamma = 2*0.01
        f = h.File(args.metacal_cat, 'r')
        index =  f['index']
        select = np.array(index['select'])
        select_1p = np.array(index['select_1p']); select_1m = np.array(index['select_1m'])
        select_2p = np.array(index['select_2p']); select_2m = np.array(index['select_2m']) 

        names=['JKR', 'ANGBIN','THETA', 'TAUO0', 'VAR_TAUO0', 'TAUO2', 'VAR_TAUO2','TAUO5', 'VAR_TAUO5']
        forms = ['i4', 'i4', 'f8',  'f8',  'f8', 'f8',  'f8',  'f8', 'f8']
        dtype = dict(names = names, formats=forms)
        nrows = 20
        outdata = np.recarray((nrows, ), dtype=dtype)

        n = h.File(args.nz_source, 'r')
        zbin_array = np.array(n['nofz/zbin'])
        nbins = 4
        for bin_c in range(nbins):
            print('Starting bin!',  bin_c)
            ind = np.where( zbin_array==bin_c )[0]
            ind_1p = np.where(np.array(n['nofz/zbin_1p'])==bin_c)
            ind_1m = np.where(np.array(n['nofz/zbin_1m'])==bin_c)
            ind_2p = np.where(np.array(n['nofz/zbin_2p'])==bin_c)
            ind_2m = np.where(np.array(n['nofz/zbin_2m'])==bin_c)

            njk = 100
            #jkindexes_stars = jk_kmeans(data_stars['ra'], data_stars['dec'] ,njk,plot=False)
            jkindexes_gals = jk_kmeans(data_gal['ra'][select][ind], data_gal['dec'][select][ind] ,njk,plot=False)
    
            for jkidx in range(njk):
                print("running jackkniffe region",  jkidx)
                booljk = [jkindexes_gals!=jkidx ] 
                R11s = (data_gal['e_1'][select_1p][ind_1p][booljk].mean() -
                        data_gal['e_1'][select_1m][ind_1m][booljk].mean() )/dgamma
                R22s = (data_gal['e_2'][select_2p][ind_2p][booljk].mean() -
                        data_gal['e_2'][select_2m][ind_2m][booljk].mean() )/dgamma
                Rs = [R11s, R22s] 
                #tau0, tau2, tau5= measure_tau(data_stars[jkindexes_stars!=jkidx ], data_gal[select][ind][booljk], Rs,   mod=args.mod, obs=args.obs)
                tau0, tau2, tau5= measure_tau(data_stars[jkindexes_stars!=jkidx ], data_gal[select][ind][booljk], Rs,   mod=args.mod, obs=args.obs)
                jkrarr =  np.array([jkidx]*nrows)
                angarr = np.arange(nrows)
                thetaarr = np.exp(tau0.meanlogr)
                if args.xim:
                    tau0arr = tau0.xim; tau2arr = tau2.xim;  tau5arr = tau5.xim;
                else:
                    tau0arr = tau0.xip; tau2arr = tau2.xip;  tau5arr = tau5.xip;
                    vartau0arr = 2*tau0.varxi; vartau2arr = 2*tau2.varxi; vartau5arr = 2*tau5.varxi;
                    array_list = [jkrarr, angar, thetaarr, tau0arr,
                                  vartau0arr, tau2arr, vartau2arr, tau5arr,
                                  vartau5arr, ]
                for array, name in zip(array_list, names): outdata[name] = array
                write_fit(outdata, names, outpath+'alltausp_jk_'+str(bin_c+1)+'_'+str(bin_c+1)+'.fits' )


            
    galkeys = ['ra','dec','e_1','e_2','R11','R22']
    data_galaxies =  read_h5(args.metacal_cat, 'catalog/metacal/unsheared',  galkeys )
    print("Total objects in catalog:", len(data_galaxies))
    dgamma = 2*0.01
    f = h.File(args.metacal_cat, 'r')
    index =  f['index']
    select = np.array(index['select'])
    select_1p = np.array(index['select_1p'])
    select_1m = np.array(index['select_1m'])
    select_2p = np.array(index['select_2p']) #added by Lucas: 
    select_2m = np.array(index['select_2m']) #added by Lucas

    R11s = (data_galaxies['e_1'][select_1p].mean() - data_galaxies['e_1'][select_1m].mean() )/dgamma
    R22s = (data_galaxies['e_2'][select_2p].mean() - data_galaxies['e_2'][select_2m].mean() )/dgamma
    Rs = [R11s, R22s]
    data_galaxies = data_galaxies[select]
    del f, index, select_1p, select_1m, select_2p,  select_2m, select

  
    names=['JKR', 'ANGBIN','THETA', 'TAUO0', 'VAR_TAUO0', 'TAUO2', 'VAR_TAUO2','TAUO5', 'VAR_TAUO5']
    forms = ['i4', 'i4', 'f8',  'f8',  'f8', 'f8',  'f8',  'f8', 'f8']
    dtype = dict(names = names, formats=forms)
    nrows = 20
    outdata = np.recarray((nrows, ), dtype=dtype)

    njk = 50
    #jkindexes_stars = jk_kmeans(data_stars['ra'], data_stars['dec'] ,njk,plot=False)
    jkindexes_gals = jk_kmeans(data_galaxies['ra'], data_galaxies['dec'] ,njk,plot=False)
    
    for jkidx in range(njk):
        print("running jackkniffe region",  jkidx)
        booljk = [jkindexes_gals!=jkidx ]  
        #tau0, tau2, tau5= measure_tau(data_stars[jkindexes_stars!=jkidx ], data_galaxies[booljk], Rs,   mod=args.mod, obs=args.obs)
        tau0, tau2, tau5= measure_tau(data_stars, data_galaxies[booljk], Rs,   mod=args.mod, obs=args.obs)
        jkrarr =  np.array([jkidx]*nrows)
        angarr = np.arange(nrows)
        thetaarr = np.exp(tau0.meanlogr)
        if args.xim:
            tau0arr = tau0.xim; tau2arr = tau2.xim;  tau5arr = tau5.xim;
        else:
            tau0arr = tau0.xip; tau2arr = tau2.xip;  tau5arr = tau5.xip;
        vartau0arr = 2*tau0.varxi; vartau2arr = 2*tau2.varxi; vartau5arr = 2*tau5.varxi;
        array_list = [jkrarr, angar, thetaarr, tau0arr, vartau0arr,
                      tau2arr, vartau2arr, tau5arr, vartau5arr, ]
        for array, name in zip(array_list, names): outdata[name] = array
        write_fit(outdata, names, outpath + args.filename)
    
if __name__ == "__main__":
    main()
