def measure_rho(data, max_sep, tag=None, use_xy=False, alt_tt=False, prefix='piff'):
    """Compute the rho statistics
    """
    import treecorr

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
    for (cat1, cat2) in [ (decat, decat),
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
