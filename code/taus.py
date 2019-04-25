import os
today = '22-04-19_'

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Produce Tau correlations, i.e correlation among galaxies and reserved stars')
    
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
    args = parser.parse_args()
    return args
    
def measure_tau(data_stars, data_galaxies, max_sep=250, sep_units='arcmin', prefix='piff', mod=True):
    """Compute the tau statistics
    """
    import treecorr
    import numpy as np
    import gc
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
        e1gal = e1gal - np.array(np.mean(e1gal))
        e2gal = e2gal - np.array(np.mean(e2gal))

        
    ra = data_stars['ra']
    dec = data_stars['dec']
    print('ra = ',ra)
    print('dec = ',dec)
    ragal = data_galaxies['ra']
    decgal = data_galaxies['dec']
    print('ragal = ',ragal)
    print('decgal = ',decgal)

    del data_stars, data_galaxies
    gc.collect()
    
    
    ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=p_e1, g2=p_e2)
    decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
    wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1, g2=w2)
    egal_cat = treecorr.Catalog(ra=ragal, dec=decgal, ra_units='deg', dec_units='deg', g1=e1gal, g2=e2gal)
    ecat.name = 'ecat'
    decat.name = 'decat'
    wcat.name = 'wcat'
    egal_cat.name = 'egal_cat'
    
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

def main():
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    #sys.path.insert(0, '/global/cscratch1/sd/alsina/alpha-beta-gamma/code/src')
    
    import numpy as np
    from read_psf_cats import read_data, toList, read_metacal
    from astropy.io import fits
    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    ##Format of the fit file output
    names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
    forms = ['i4', 'i4', 'i4',  'f8',  'f8']
    dtype = dict(names = names, formats=forms)
    nrows = 20
    outdata = np.recarray((nrows, ), dtype=dtype)
    namesout=['TAU0P', 'TAU2P', 'TAU5P', 'TAU0M','TAU2M', 'TAU5M']

  
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
        print('Starting Tomography!')
        galkeys = ['ra','dec','e_1','e_2','R11','R22']
        nbins = 4
    
        for bin_c in range(nbins):
            print('Starting bin!',  bin_c)
            data_gal = read_metacal(args.metacal_cat,  galkeys,  zbin=bin_c,  nz_source_file=args.nz_source)
            tau0, tau2, tau5= measure_tau(data_stars, data_gal,  mod=args.mod)
            tau0parr = tau0.xip; tau2parr = tau2.xip;  tau5parr = tau5.xip;
            tau0marr = tau0.xim; tau2marr = tau2.xim;  tau5marr = tau5.xim;
            vartau0arr = 2*tau0.varxi; vartau2arr = 2*tau2.varxi; vartau5arr = 2*tau5.varxi;
            taus = [tau0parr, tau2parr, tau5parr, tau0marr, tau2marr,
                    tau5marr]
            vares = [vartau0arr, vartau2arr, vartau5arr, vartau0arr,
                     vartau2arr, vartau5arr]
            
            for i, nam in enumerate(namesout):
                covmat = np.diag(vares[i])
                hdu = fits.PrimaryHDU()
                hdul = fits.HDUList([hdu])
                covmathdu = fits.ImageHDU(covmat, name='COVMAT')
                hdul.insert(1, covmathdu)
          
                angarray = np.exp(tau0.meanlogr)
                valuearray =  np.array(taus[i])
                bin1array = np.array([bin_c]*nrows)
                bin2array = np.array([bin_c]*nrows)
                angbinarray = np.arange(nrows)
                array_list = [bin1array, bin2array, angbinarray, valuearray,  angarray ]
                for array, name in zip(array_list, names): outdata[name] = array 
                corrhdu = fits.BinTableHDU(outdata, name=nam)
                hdul.insert(2, corrhdu)
    
                hdul.writeto(outpath + nam +'_bin_' + str(bin_c) +  '.fits', clobber=True)
                

            
    galkeys = ['ra','dec','e_1','e_2','R11','R22']
    data_galaxies =  read_metacal(args.metacal_cat,  galkeys )
    print("Total objects in catalog:", len(data_galaxies))
    
    tau0, tau2, tau5= measure_tau(data_stars, data_galaxies, mod=args.mod)
    angarr = np.arange(nrows)
    thetaarr = np.exp(tau0.meanlogr)
    tau0marr = tau0.xim; tau2marr = tau2.xim;  tau5marr = tau5.xim;
    tau0parr = tau0.xip; tau2parr = tau2.xip;  tau5parr = tau5.xip;
    vartau0arr = 2*tau0.varxi; vartau2arr = 2*tau2.varxi; vartau5arr = 2*tau5.varxi;

    taus = [tau0parr, tau2parr, tau5parr, tau0marr, tau2marr, tau5marr ]
    vares = [vartau0arr, vartau2arr, vartau5arr, vartau0arr, vartau2arr, vartau5arr]
    for i, nam in enumerate(namesout):
        covmat = np.diag(vares[i])
        hdu = fits.PrimaryHDU()
        hdul = fits.HDUList([hdu])
        covmathdu = fits.ImageHDU(covmat, name='COVMAT')
        hdul.insert(1, covmathdu)
        
        angarray = np.exp(tau0.meanlogr)
        valuearray =  np.array(taus[i])
        bin1array = np.array([ -999]*nrows)
        bin2array = np.array([ -999]*nrows)
        angbinarray = np.arange(nrows)
        array_list = [bin1array, bin2array, angbinarray, valuearray,  angarray ]
        for array, name in zip(array_list, names): outdata[name] = array 

        corrhdu = fits.BinTableHDU(outdata, name=nam)
        hdul.insert(2, corrhdu)
    
        hdul.writeto(outpath + nam + '.fits', clobber=True)
    
if __name__ == "__main__":
    main()
