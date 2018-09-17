def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Alpha beta gamma test')
    
    parser.add_argument('--metacal_cat',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3_master/Y3_mastercat_v2_6_20_18_subsampled.h5',
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--piff_cat',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3a1-v29',
                        help='Full Path to the Only stars Piff catalog')
    parser.add_argument('--exps_file',
                        default='/home/dfa/sobreira/alsina/DESWL/psf/ally3.grizY',
                        help='list of exposures (in lieu of separate exps)')
    parser.add_argument('--bands', default='grizY', type=str,
                         help='Limit to the given bands')
    parser.add_argument('--use_reserved', default=False,
                        action='store_const', const=True,
                        help='just use the objects with the RESERVED flag')
    parser.add_argument('--frac', default=1., type=float,
                        help='Choose a random fraction of the input stars')
    
    
    args = parser.parse_args()

    return args


def measure_rho(data, max_sep, tag=None, use_xy=False, alt_tt=False, prefix='piff', lucas=False):
    """Compute the rho statistics
    """
    import treecorr

    e1 = data['obs_e1']
    e2 = data['obs_e2']
    T = data['obs_T']
    p_e1 = data[prefix+'_e1']
    p_e2 = data[prefix+'_e2']
    p_T = data[prefix+'_T']
    m = data['mag']

    mlt20 = True
    if mlt20:
        e1 = e1[m<20]
        e2 = e2[m<20]
        T = T[m<20]
        p_e1 = p_e1[m<20]
        p_e2 = p_e2[m<20]
        p_T = p_T[m<20]

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
        if mlt20:
            x = x[m<20]
            y = y[m<20]
        print('x = ',x)
        print('y = ',y)

        ecat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec', g1=e1, g2=e2)
        decat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec', g1=de1, g2=de2)
        dtcat = treecorr.Catalog(x=x, y=y, x_units='arcsec', y_units='arcsec',
                                 k=dt, g1=dt*e1, g2=dt*e2)
    else:
        ra = data['ra']
        dec = data['dec']
        if mlt20:
            ra = ra[m<20]
            dec = dec[m<20]
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
    )

    if lucas:
        bin_config['min_sep'] = 2.5
        bin_config['max_sep'] = 250.
        bin_config['nbins'] = 20
        del bin_config['bin_size']

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



def main():
    import numpy as np
    import h5py as h
    from read_psf_cats import read_data
    from read_psf_cats import toList
    
    args = parse_args()

    f = h.File(args.metacal_cat, 'r')
    print(f.keys())
    catalog =  f['catalog']
    print('Master catalog folders', catalog.keys())
    metacal =  catalog['metacal']
    print('Metacal catalog folders',metacal.keys())
    #data = f.get('catalog')
    #dataset1 = np.array(data)
    #print(data['metacal'])
    unsheared =  metacal['unsheared']
    print('Unsheared Metacal catalog fields',unsheared.keys())

    keys = ['ra', 'dec','obs_e1', 'obs_e2', 'obs_T',
            'piff_e1', 'piff_e2', 'piff_T']
 
    exps = toList(args.exps_file)
    data, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved)

if __name__ == "__main__":
    main()
