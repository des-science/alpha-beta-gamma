import os
def write_stats(stat_file, rho0, rho1, rho2, rho3, rho4, rho5, corr_tt=None):
    import json

    stats = [
        rho0.meanlogr.tolist(),
        rho0.xip.tolist(),
        rho0.xip_im.tolist(),
        rho0.xim.tolist(),
        rho0.xim_im.tolist(),
        rho0.varxi.tolist(),
        rho1.xip.tolist(),
        rho1.xip_im.tolist(),
        rho1.xim.tolist(),
        rho1.xim_im.tolist(),
        rho1.varxi.tolist(),
        rho2.xip.tolist(),
        rho2.xip_im.tolist(),
        rho2.xim.tolist(),
        rho2.xim_im.tolist(),
        rho2.varxi.tolist(),
        rho3.xip.tolist(),
        rho3.xip_im.tolist(),
        rho3.xim.tolist(),
        rho3.xim_im.tolist(),
        rho3.varxi.tolist(),
        rho4.xip.tolist(),
        rho4.xip_im.tolist(),
        rho4.xim.tolist(),
        rho4.xim_im.tolist(),
        rho4.varxi.tolist(),
        rho5.xip.tolist(),
        rho5.xip_im.tolist(),
        rho5.xim.tolist(),
        rho5.xim_im.tolist(),
        rho5.varxi.tolist(),
    ]
    if corr_tt is not None:
        stats.extend([
            corr_tt.xi.tolist(),
            corr_tt.varxi.tolist()
        ])
    #print('stats = ',stats)
    print('stat_file = ',stat_file)
    with open(stat_file,'w') as fp:
        json.dump([stats], fp)
    print('Done writing ',stat_file)

def write_cross_stats(stat_file, sigma0, sigma2, sigma5, corr_tt=None):
    import json

    stats = [
        sigma0.meanlogr.tolist(),
        sigma0.xip.tolist(),
        sigma0.xip_im.tolist(),
        sigma0.xim.tolist(),
        sigma0.xim_im.tolist(),
        sigma0.varxi.tolist(),
        sigma2.xip.tolist(),
        sigma2.xip_im.tolist(),
        sigma2.xim.tolist(),
        sigma2.xim_im.tolist(),
        sigma2.varxi.tolist(),
        sigma5.xip.tolist(),
        sigma5.xip_im.tolist(),
        sigma5.xim.tolist(),
        sigma5.xim_im.tolist(),
        sigma5.varxi.tolist(),
    ]
    if corr_tt is not None:
        stats.extend([
            corr_tt.xi.tolist(),
            corr_tt.varxi.tolist()
        ])
    #print('stats = ',stats)
    print('stat_file = ',stat_file)
    with open(stat_file,'w') as fp:
        json.dump([stats], fp)
    print('Done writing ',stat_file)

def measure_rho(data, max_sep, tag=None, use_xy=False, alt_tt=False, prefix='piff'):
    """Compute the rho statistics
    """
    import treecorr
    import numpy as np

    e1 = data['obs_e1']
    e2 = data['obs_e2']
    T = data['obs_T']
    p_e1 = data[prefix+'_e1']
    p_e2 = data[prefix+'_e2']
    p_T = data[prefix+'_T']

    de1 = e1-p_e1
    de2 = e2-p_e2
    dt = (T-p_T)/T
    print('mean e = ',np.mean(e1),np.mean(e2))
    print('mean T = ',np.mean(T))
    print('mean de = ',np.mean(de1),np.mean(de2))
    print('mean dT = ',np.mean(T-p_T))
    print('mean dT/T = ',np.mean(dt))

    if use_xy:
        x = data['fov_x']
        y = data['fov_y']
        print('x = ',x)
        print('y = ',y)

        ecat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec', g1=e1, g2=e2)
        decat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec', g1=de1, g2=de2)
        dtcat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec',
                                 k=dt, g1=dt*e1, g2=dt*e2)
    else:
        ra = data['ra']
        dec = data['dec']
        print('ra = ',ra)
        print('dec = ',dec)

        ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=e1, g2=e2)
        decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
        dtcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg',
                                 k=dt, g1=dt*e1, g2=dt*e2)
    ecat.name = 'ecat'
    decat.name = 'decat'
    dtcat.name = 'dtcat'
    if tag is not None:
        for cat in [ ecat, decat, dtcat ]:
            cat.name = tag + ":"  + cat.name

    bin_config = dict(
        sep_units = 'arcmin',
        bin_slop = 0.1,

        min_sep = 0.5,
        max_sep = max_sep,
        bin_size = 0.2,

        #min_sep = 2.5,
        #max_sep = 250,
        #nbins = 20,
    )

    results = []
    for (cat1, cat2) in [(ecat, ecat), 
                         (decat, decat),
                          (ecat, decat),
                          (dtcat, dtcat),
                          (decat, dtcat),
                          (ecat, dtcat) ]:
        print('Doing correlation of %s vs %s'%(cat1.name, cat2.name))

        rho = treecorr.GGCorrelation(bin_config, verbose=2)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        print('mean xi+ = ',rho.xip.mean())
        print('mean xi- = ',rho.xim.mean())
        results.append(rho)

    if alt_tt:
        print('Doing alt correlation of %s vs %s'%(dtcat.name, dtcat.name))

        rho = treecorr.KKCorrelation(bin_config, verbose=2)
        rho.process(dtcat)
        results.append(rho)

    return results

def measure_cross_rho(data_stars, data_galaxies, max_sep, tag=None, use_xy=False, alt_tt=False, prefix='piff'):
    """Compute the rho statistics
    """
    import treecorr
    import numpy as np

    e1 = data_stars['obs_e1']
    e2 = data_stars['obs_e2']
    T = data_stars['obs_T']
    p_e1 = data_stars[prefix+'_e1']
    p_e2 = data_stars[prefix+'_e2']
    p_T = data_stars[prefix+'_T']

    e1gal =  data_galaxies['e_1']
    e2gal =  data_galaxies['e_2']

    de1 = e1-p_e1
    de2 = e2-p_e2
    dt = (T-p_T)/T
    print('mean e = ',np.mean(e1),np.mean(e2))
    print('mean T = ',np.mean(T))
    print('mean de = ',np.mean(de1),np.mean(de2))
    print('mean dT = ',np.mean(T-p_T))
    print('mean dT/T = ',np.mean(dt))

    ra = data_stars['ra']
    dec = data_stars['dec']
    print('ra = ',ra)
    print('dec = ',dec)
    ragal = data_galaxies['ra']
    decgal = data_galaxies['dec']
    print('ragal = ',ragal)
    print('decgal = ',decgal)
    
    ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=e1, g2=e2)
    decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
    dtcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg',
                             k=dt, g1=dt*e1, g2=dt*e2)
    egal_cat = treecorr.Catalog(ra=ragal, dec=decgal, ra_units='deg', dec_units='deg', g1=e1gal, g2=e2gal)

    ecat.name = 'ecat'
    decat.name = 'decat'
    dtcat.name = 'dtcat'
    egal_cat.name = 'egal_cat'
    
    if tag is not None:
        for cat in [ ecat, decat, dtcat ]:
            cat.name = tag + ":"  + cat.name

    bin_config = dict(
        sep_units = 'arcmin',
        bin_slop = 0.1,

        min_sep = 0.5,
        max_sep = max_sep,
        bin_size = 0.2,

        #min_sep = 2.5,
        #max_sep = 250,
        #nbins = 20,
    )

    results = []
    for (cat1, cat2) in [(egal_cat, ecat), 
                         (egal_cat, decat),
                          (egal_cat, dtcat) ]:
        print('Doing correlation of %s vs %s'%(cat1.name, cat2.name))

        rho = treecorr.GGCorrelation(bin_config, verbose=2)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        print('mean xi+ = ',rho.xip.mean())
        print('mean xi- = ',rho.xim.mean())
        results.append(rho)

    if alt_tt:
        print('Doing alt correlation of %s vs %s'%(dtcat.name, dtcat.name))

        rho = treecorr.KKCorrelation(bin_config, verbose=2)
        rho.process(dtcat)
        results.append(rho)

    return results

def band_combinations(bands, single=True, combo=True,  allcombo=True):
    if(allcombo):
        if single:
            use_bands = [ [b] for b in bands ]
        else:
            use_bands = []

        if combo:
            if 'r' in bands and 'i' in bands:
                use_bands.append(['r', 'i'])
            if 'r' in bands and 'i' in bands and 'z' in bands:
                use_bands.append(['r', 'i', 'z'])
            if 'g' in bands and 'r' in bands and 'i' in bands and 'z' in bands:
                use_bands.append(['g', 'r', 'i', 'z'])
    else:
        letters = [k for k in bands]
        use_bands = []
        use_bands.append(letters)

    print('use_bands = ',use_bands)
    print('tags = ',[ ''.join(band) for band in use_bands ])
    return use_bands
def do_canonical_stats(data, bands, tilings, outpath, prefix='piff', name='all', alt_tt=False, bandcombo=True):
    import numpy as np 
    print('Start CANONICAL: ',prefix,name)
    # Measure the canonical rho stats using all pairs:
    use_bands = band_combinations(bands, allcombo=bandcombo)
    for band in use_bands:
        print('band ',band)
        mask = np.in1d(data['band'],band)
        print('sum(mask) = ',np.sum(mask))
        print('len(data[mask]) = ',len(data[mask]))
        tag = ''.join(band)
        stats = measure_rho(data[mask], max_sep=300, tag=tag, prefix=prefix, alt_tt=alt_tt)
        stat_file = os.path.join(outpath, "rho_%s_%s.json"%(name,tag))
        write_stats(stat_file,*stats)

def do_cross_stats(data_stars, data_galaxies,  bands, tilings, outpath, prefix='piff', name='all', alt_tt=False, bandcombo=True):
    import numpy as np 
    print('Start CANONICAL: ',prefix,name)
    # Measure the canonical rho stats using all pairs:
    use_bands = band_combinations(bands, allcombo=bandcombo)
    for band in use_bands:
        print('band ',band)
        mask_stars = np.in1d(data_stars['band'],band)
        #mask_galaxies = np.in1d(data_galaxies['band'],band)
        print('sum(mask) = ',np.sum(mask_stars))
        print('len(data[mask]) = ',len(data_stars[mask_stars]))
        tag = ''.join(band)
        #stats = measure_cross_rho(data_stars[mask_stars], data_galaxies[mask_galaxies], max_sep=300, tag=tag, prefix=prefix, alt_tt=alt_tt)
        stats = measure_cross_rho(data_stars[mask_stars], data_galaxies, max_sep=300, tag=tag, prefix=prefix, alt_tt=alt_tt)
        stat_file = os.path.join(outpath, "sigma_%s_%s.json"%(name,tag))
        write_cross_stats(stat_file,*stats)
