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

def measure_rho(data, max_sep, tag=None, use_xy=False, alt_tt=False, prefix='piff', mod=True):
    """Compute the rho statistics
    """
    import treecorr
    import numpy as np

    e1 = data['obs_e1']
    e2 = data['obs_e2']
    
    p_e1 = data[prefix+'_e1']
    p_e2 = data[prefix+'_e2']
    
    de1 = e1-p_e1
    de2 = e2-p_e2

    
    T = data['obs_T']
    p_T = data[prefix+'_T']
    dt = (T-p_T)/T
    w1 = p_e1*dt
    w2 = p_e2*dt
    
    #Modified ellipticities
    if(mod):
        e1 = e1 - np.array(np.mean(e1))
        e2 = e2 - np.array(np.mean(e2))
        de1 = de1 - np.array(np.mean(de1))
        de2 = de2 - np.array(np.mean(de2))
        w1 = w1 - np.array(np.mean(w1))
        w2 = w2 - np.array(np.mean(w2))
    
    if use_xy:
        x = data['fov_x']
        y = data['fov_y']
        print('x = ',x)
        print('y = ',y)

        ecat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec', g1=p_e1, g2=p_e2)
        decat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec', g1=de1, g2=de2)
        wcat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec',
                                 k=dt, g1=w1, g2=w2)
    else:
        ra = data['ra']
        dec = data['dec']
        print('ra = ',ra)
        print('dec = ',dec)

        ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=p_e1, g2=p_e2)
        decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
        wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg',
                                 k=dt, g1=w1, g2=w2)
    ecat.name = 'ecat'
    decat.name = 'decat'
    wcat.name = 'wcat'
    if tag is not None:
        for cat in [ ecat, decat, wcat ]:
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

    if alt_tt:
        print('Doing alt correlation of %s vs %s'%(wcat.name, wcat.name))

        rho = treecorr.KKCorrelation(bin_config, verbose=2)
        rho.process(wcat)
        results.append(rho)

    return results

def measure_cross_rho(data_stars, data_galaxies, max_sep, tag=None, use_xy=False, alt_tt=False, prefix='piff', mod=True):
    """Compute the rho statistics
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

    w1 = p_e1*dt
    w2 = p_e2*dt

    #Modified ellipticities reserved stars
    e1 = e1 - np.array(np.mean(e1))
    e2 = e2 - np.array(np.mean(e2))
    de1 = de1 - np.array(np.mean(de1))
    de2 = de2 - np.array(np.mean(de2))
    w1 = w1 - np.array(np.mean(w1))
    w2 = w2 - np.array(np.mean(w2))

    e1gal = data_galaxies['e_1']
    e2gal = data_galaxies['e_2']
    R11 =  data_galaxies['R11']
    R22 =  data_galaxies['R22']

    #Modified ellipticities galaxies
    if (mod):
        e1gal = (e1gal - np.array(np.mean(e1gal)))/(np.mean(R11)) 
        e2gal = (e2gal - np.array(np.mean(e2gal)))/(np.mean(R22)) 

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
    wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg',
                             k=dt, g1=w1, g2=w2)
    egal_cat = treecorr.Catalog(ra=ragal, dec=decgal, ra_units='deg', dec_units='deg', g1=e1gal, g2=e2gal)

    ecat.name = 'ecat'
    decat.name = 'decat'
    wcat.name = 'wcat'
    egal_cat.name = 'egal_cat'
    
    if tag is not None:
        for cat in [ ecat, decat, wcat ]:
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

    if alt_tt:
        print('Doing alt correlation of %s vs %s'%(wcat.name, wcat.name))

        rho = treecorr.KKCorrelation(bin_config, verbose=2)
        rho.process(wcat)
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
def do_canonical_stats(data, bands, tilings, outpath, prefix='piff', name='all', alt_tt=False, bandcombo=True, mod=True):
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
        stats = measure_rho(data[mask], max_sep=300, tag=tag, prefix=prefix, alt_tt=alt_tt, mod=mod)
        stat_file = os.path.join(outpath, "rho_%s_%s.json"%(name,tag))
        write_stats(stat_file,*stats)

def do_cross_stats(data_stars, data_galaxies,  bands, tilings, outpath, prefix='piff', name='all', alt_tt=False, bandcombo=True, mod=True):
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
        stats = measure_cross_rho(data_stars[mask_stars], data_galaxies, max_sep=300, tag=tag, prefix=prefix, alt_tt=alt_tt,  mod=mod)
        stat_file = os.path.join(outpath, "tau_%s_%s.json"%(name,tag))
        write_cross_stats(stat_file,*stats)
