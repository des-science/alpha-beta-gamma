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
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/',
                        help='location of the output of the files')
    parser.add_argument('--tomo', default=False,
                        action='store_const', const=True,
                        help='Run all tomographic correlations')
    parser.add_argument('--nz_source',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3_master/nz_source_zbin.h5',
                        help='Full Path to the Only stars Piff catalog')
    parser.add_argument('--filename', default='alltaus_4jk.fits', help='Name of the fit file where info of dxip will be saved ')
    
    
    args = parser.parse_args()

    return args

def jk_kmeans(ra_sam, dec_sam, ra,dec,njk,plot=False):
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
        from astropy.coordinates import SkyCoord, Angle
        from astropy import units
        radec = np.zeros((len(ra),2)); radec_sam = np.zeros((len(ra_sam),2))
        radec[:,0] = ra; radec_sam[:,0] = ra_sam
        radec[:,1] = dec; radec_sam[:,1] = dec_sam
        km = kmeans_radec.kmeans_sample(radec_sam,njk,maxiter=500,tol=1e-05)
        jk = km.find_nearest(radec)
        if not km.converged:
                print 'k means did not converge'
        if plot:
            coords = SkyCoord(ra=ra, dec=dec, unit='degree')
            ra = coords.ra.wrap_at(180 * units.deg)
            dec = coords.dec
            plt.figure()
            plt.scatter(ra,dec,c=jk,lw=0,cmap='Paired',rasterized=True)
            plt.xlabel(r'RA',fontsize=12)
            plt.ylabel(r'Dec',fontsize=12)
            plt.tight_layout()
            plt.savefig('jk_kmeans.png')
        return jk

def measure_tau(data_stars, data_galaxies, max_sep=250, sep_units='arcmin', prefix='piff', mod=True):
    """Compute the tau statistics
    """
    import treecorr
    import numpy as np
    e1 = data_stars['obs_e1']
    e2 = data_stars['obs_e2']
    p_e1 = data_stars[prefix+'_e1']
    p_e2 = data_stars[prefix+'_e2']
    T = data_stars['obs_T']
    p_T = data_stars[prefix+'_T']

    de1 = e1-p_e1
    de2 = e2-p_e2
    dt = (T-p_T)/T

    #w1 = p_e1*dt
    #w2 = p_e2*dt
    w1 = e1*dt
    w2 = e2*dt

    e1gal = data_galaxies['e_1']
    e2gal = data_galaxies['e_2']
    
    #Modified ellipticities reserved stars and galaxies
    if(mod):
        p_e1 = p_e1 - np.array(np.mean(p_e1))
        p_e2 = p_e2 - np.array(np.mean(p_e2))
        de1 = de1 - np.array(np.mean(de1))
        de2 = de2 - np.array(np.mean(de2))
        w1 = w1 - np.array(np.mean(w1))
        w2 = w2 - np.array(np.mean(w2))
        e1gal = (e1gal - np.array(np.mean(e1gal)))
        e2gal = (e2gal - np.array(np.mean(e2gal)))

        
    ra = data_stars['ra']
    dec = data_stars['dec']
    print('ra = ',ra)
    print('dec = ',dec)
    ragal = data_galaxies['ra']
    decgal = data_galaxies['dec']
    print('ragal = ',ragal)
    print('decgal = ',decgal)
    
    
    ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=p_e1, g2=p_e2)
    decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
    wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1, g2=w2)
    egal_cat = treecorr.Catalog(ra=ragal, dec=decgal, ra_units='deg', dec_units='deg', g1=e1gal, g2=e2gal)

    bin_config = dict( sep_units = sep_units, min_sep = 2.5, max_sep = 250, nbins = 20,)
    #bin_config = dict(sep_units = 'degrees', bin_slop = 0.1, min_sep = 0.5, max_sep = max_sep, bin_size = 0.2)
  
    results = []
    for (cat1, cat2) in [(egal_cat, ecat), 
                         (egal_cat, decat),
                          (egal_cat, wcat) ]:
        print('Doing correlation of %s vs %s'%(cat1.name, cat2.name))

        rho = treecorr.GGCorrelation(bin_config, verbose=2)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        print('mean xi+ = ',rho.xip.mean())
        print('mean xi- = ',rho.xim.mean())
        results.append(rho)
        
    return results

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
    from read_psf_cats import read_data, toList, read_metacal
    import h5py as h
    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
        
    names=['JKR','ANGBIN','THETA','TAU0P','TAU0M','VAR_TAU0','TAU2P','TAU2M','VAR_TAU2','TAU5P','TAU5M', 'VAR_TAU5']
    forms = ['i4', 'i4', 'f8',  'f8',  'f8', 'f8',  'f8',  'f8', 'f8',  'f8',  'f8', 'f8']
    dtype = dict(names = names, formats=forms)
    nrows = 20
    outdata = np.recarray((nrows, ), dtype=dtype)
    
  
    #Reading Mike stars catalog
    keys = ['ra', 'dec','obs_e1', 'obs_e2', 'obs_T',
            'piff_e1', 'piff_e2', 'piff_T', 'mag']
 
    exps = toList(args.exps_file)
    data_sam, bands, tilings = read_data(exps, args.piff_cat , keys,
                                          limit_bands=args.bands,
                                          use_reserved=args.use_reserved, frac=0.01)
    data_stars, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved)
    print("Objects",  len(data_stars))
    data_stars = data_stars[data_stars['mag']<20]
    print("Objects with magnitude <20",  len(data_stars))
    
    
    if(args.tomo):
        print('Starting Tomography!')
        galkeys = ['ra','dec','e_1','e_2','R11','R22']
        nbins = 4
        for bin_c in range(nbins):
            print('Starting bin!',  bin_c)
            data_gal = read_metacal(args.metacal_cat,  galkeys,  zbin=bin_c,  nz_source_file=args.nz_source)
            
            njk = 4
            ##TODO generate km first an later finnearest,
            jkindexes_gals = jk_kmeans(data_sam['ra'], data_sam['dec'], data_gal['ra'], data_gal['dec'],njk)
    
            for jkidx in range(njk):
                print("running jackkniffe region",  jkidx)
                booljk = [jkindexes_gals!=jkidx ] 
                tau0, tau2, tau5= measure_tau(data_stars, data_gal[booljk], mod=args.mod)
                jkrarr =  np.array([jkidx]*nrows)
                angarr = np.arange(nrows)
                thetaarr = np.exp(tau0.meanlogr)
                tau0parr = tau0.xip; tau2parr = tau2.xip;  tau5parr = tau5.xip;
                tau0marr = tau0.xim; tau2marr = tau2.xim;  tau5marr = tau5.xim;
                vartau0arr = 2*tau0.varxi; vartau2arr = 2*tau2.varxi; vartau5arr = 2*tau5.varxi;
                array_list = [jkrarr, angar, thetaarr,
                              tau0parr,tau0marr, vartau0arr,
                              tau2parr,tau2marr, vartau2arr,
                              tau5parr,tau5marr, vartau5arr, ]
                for array, name in zip(array_list, names): outdata[name] = array
                write_fit(outdata, names, outpath+'alltaus_4jk_'+str(bin_c+1)+'_'+str(bin_c+1)+'.fits' )


            
    galkeys = ['ra','dec','e_1','e_2','R11','R22']
    data_galaxies =  read_metacal(args.metacal_cat,  galkeys )
    print("Total objects in catalog:", len(data_galaxies))
 
   
    njk = 4
 
    jkindexes_gals = jk_kmeans(data_sam['ra'], data_sam['dec'], data_galaxies['ra'], data_galaxies['dec'],njk)
    
    for jkidx in range(njk):
        print("running jackkniffe region",  jkidx)
        booljk = [jkindexes_gals!=jkidx ]  
        tau0, tau2, tau5= measure_tau(data_stars, data_galaxies[booljk], mod=args.mod)
        jkrarr =  np.array([jkidx]*nrows)
        angarr = np.arange(nrows)
        thetaarr = np.exp(tau0.meanlogr)
        tau0parr = tau0.xip; tau2parr = tau2.xip;  tau5parr = tau5.xip;
        tau0marr = tau0.xim; tau2marr = tau2.xim;  tau5marr = tau5.xim;
        vartau0arr = 2*tau0.varxi; vartau2arr = 2*tau2.varxi; vartau5arr = 2*tau5.varxi;
        array_list = [jkrarr, angarr, thetaarr, tau0parr,tau0marr,
                      vartau0arr, tau2parr,tau2marr, vartau2arr,
                      tau5parr,tau5marr, vartau5arr, ]
        for array, name in zip(array_list, names): outdata[name] = array
        write_fit(outdata, names, outpath + args.filename)
    
if __name__ == "__main__":
    main()
