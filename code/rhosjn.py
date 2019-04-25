import numpy as np
import os
import kmeans_radec
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

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
    parser.add_argument('--obs', default=False,
                        action='store_const', const=True,
                        help='Use e_obs instead of e_piff to calculate modified rho stats')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations',
                        help='location of the output of the files')
    parser.add_argument('--filename', default='allrhos_4jk.fits', help='Name of the fit file where info of dxip will be saved ')
    
    
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
def measure_rho(data, max_sep=300, sep_units='arcmin',  tag=None, prefix='piff', mod=True,  obs=False ):
    """Compute the rho statistics
    """
    import treecorr
    import numpy as np

    e1 = data['obs_e1']
    e2 = data['obs_e2']
    p_e1 = data[prefix+'_e1']
    p_e2 = data[prefix+'_e2']
    T = data['obs_T']
    p_T = data[prefix+'_T']

    de1 = e1-p_e1
    de2 = e2-p_e2
    dt = (T-p_T)/T
    w1 = p_e1*dt
    w2 = p_e2*dt
    w1obs = e1*dt 
    w2obs = e2*dt 
        
    
    #Modified ellipticities
    if(mod):
        e1 = e1 - np.array(np.mean(e1))
        e2 = e2 - np.array(np.mean(e2))       
        p_e1 = p_e1 - np.array(np.mean(p_e1))
        p_e2 = p_e2 - np.array(np.mean(p_e2))
        de1 = de1 - np.array(np.mean(de1))
        de2 = de2 - np.array(np.mean(de2))
        w1 = w1 - np.array(np.mean(w1))
        w2 = w2 - np.array(np.mean(w2))
        w1obs = w1obs - np.array(np.mean(w1obs))
        w2obs = w2obs - np.array(np.mean(w2obs))
        
    ra = data['ra']
    dec = data['dec']
    print('ra = ',ra)
    print('dec = ',dec)
    if(obs):
        ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=e1, g2=e2)
        decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
        wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1obs, g2=w2obs)
    else:
        ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=p_e1, g2=p_e2)
        decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
        #wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1, g2=w2)
        wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1obs, g2=w2obs)
    ecat.name = 'ecat'
    decat.name = 'decat'
    wcat.name = 'wcat'
    if tag is not None:
        for cat in [ ecat, decat, wcat ]:
            cat.name = tag + ":"  + cat.name

    bin_config = dict(
        sep_units = sep_units,
        nbins = 20,
        min_sep = 2.5,
        max_sep = 250,)
        
    #sep_units = 'degrees',
    '''
    bin_slop = 0.1,
    min_sep = 0.5,
    max_sep = max_sep,
    bin_size = 0.2,
    '''   
    

    results = []
    for (cat1, cat2) in [(ecat, ecat), 
                         (decat, decat),
                          (decat, ecat),
                          (wcat, wcat),
                          (decat, wcat),
                          (ecat, wcat) ]:
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
 
def main():
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    #sys.path.insert(0, '/global/cscratch1/sd/alsina/alpha-beta-gamma/code/src')

    import numpy as np
    from read_psf_cats import read_data,  toList
    
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
    data_sam, bands, tilings = read_data(exps, args.piff_cat , keys,
                                          limit_bands=args.bands,
                                          use_reserved=args.use_reserved, frac=0.01)
    data, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved)
    print("Objects",  len(data))
    data = data[data['mag']<20]
    print("Objects with magnitude <20",  len(data))
    
    names = ['JKR', 'ANGBIN','THETA', 'RHO0P', 'RHO0M', 'VAR_RHO0', 'RHO1P', 'RHO1M',
             'VAR_RHO1', 'RHO2P', 'RHO2M','VAR_RHO2', 'RHO3P', 'RHO3M','VAR_RHO3','RHO4P', 'RHO4M',
             'VAR_RHO4', 'RHO5P', 'RHO5M','VAR_RHO5']
    forms = ['i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
             'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
             'f8', 'f8', 'f8']
    dtype = dict(names = names, formats=forms)
    nrows = 20
    outdata = np.recarray((nrows, ), dtype=dtype)

    njk = 4
    jkindexes = jk_kmeans(data_sam['ra'], data_sam['dec'], data['ra'], data['dec'] ,njk,plot=True)
    #print (jkindexes)
    
    for jkidx in range(njk):
        print("running jackkniffe region",  jkidx)
        rho0, rho1, rho2, rho3, rho4, rho5 = measure_rho(data[ jkindexes!=jkidx ],  mod=args.mod, obs=args.obs)
        jkrarr =  np.array([jkidx]*nrows)
        angarr = np.arange(nrows)
        thetaarr = np.exp(rho0.meanlogr)

        rho0marr = rho0.xim; rho1marr = rho1.xim; rho2marr = rho2.xim; rho3marr = rho3.xim; rho4marr = rho4.xim; rho5marr = rho5.xim
        rho0parr = rho0.xip; rho1parr = rho1.xip; rho2parr = rho2.xip; rho3parr = rho3.xip; rho4parr = rho4.xip; rho5parr = rho5.xip
        varrho0arr = 2*rho0.varxi; varrho1arr = 2*rho1.varxi; varrho2arr = 2*rho2.varxi;
        varrho3arr = 2*rho3.varxi; varrho4arr = 2*rho4.varxi; varrho5arr = 2*rho5.varxi;
        array_list = [jkrarr, angarr, thetaarr, rho0parr, rho0marr, varrho0arr,
                      rho1parr, rho1marr, varrho1arr,rho2parr, rho2marr, varrho2arr,rho3parr, rho3marr,
                      varrho3arr,rho4parr, rho4marr, varrho4arr,rho5parr, rho5marr,
                      varrho5arr]
        for array, name in zip(array_list, names): outdata[name] = array
        write_fit(outdata, names, outpath + args.filename)
    
if __name__ == "__main__":
    main()
