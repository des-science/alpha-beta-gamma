#runing taus by quadrant.
import os
today = '08-04-19_'
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation galaxies and reserved stars')
    
    parser.add_argument('--metacal_cat',
                        #default='/home2/dfa/sobreira/alsina/catalogs/y3_master/Y3_mastercat_v2_6_20_18_subsampled.h5',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3_master/Y3_mastercat_v2_6_20_18.h5', 
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
    parser.add_argument('--outpath', default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma',
                        help='location of the output of the files')
    parser.add_argument('--tomo', default=False,
                        action='store_const', const=True,
                        help='Run all tomographic correlations')
    parser.add_argument('--nz_source',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3_master/nz_source_zbin.h5',
                        help='Full Path to the Only stars Piff catalog')
    
    
    args = parser.parse_args()

    return args

        
def main():
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    #sys.path.insert(0, '/global/cscratch1/sd/alsina/alpha-beta-gamma/code/src')
    
    import numpy as np
    from read_psf_cats import read_data, toList, read_h5
    from run_rho import  do_tau_stats
    import h5py as h
    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
        

    #Reading metacal catalog
    #galkeys = ['ra']
    #blabla =  read_h5(args.metacal_cat, 'catalog/metacal/sheared_1m',  galkeys )
        
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
    meanra = np.mean(data_stars['ra'])
    meandec = np.mean(data_stars['dec']) 
    
     
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
            R11s = (data_gal['e_1'][select_1p][ind_1p].mean() -
                data_gal['e_1'][select_1m][ind_1m].mean() )/dgamma
            R22s = (data_gal['e_2'][select_2p][ind_2p].mean() -
                  data_gal['e_2'][select_2m][ind_2m].mean() )/dgamma
            Rs = [R11s, R22s]

            patchstars = [];patchgal = []
            patchstars.append((data_stars['ra']>meanra)&(data_stars['dec']>meandec))
            patchstars.append((data_stars['ra']<meanra)&(data_stars['dec']>meandec))
            patchstars.append((data_stars['ra']<meanra)&(data_stars['dec']<meandec))
            patchstars.append((data_stars['ra']>meanra)&(data_stars['dec']<meandec))
            patchgal.append((data_gal[select][ind]['ra']>meanra)&(data_gal[select][ind]['dec']>meandec))
            patchgal.append((data_gal[select][ind]['ra']<meanra)&(data_gal[select][ind]['dec']>meandec))
            patchgal.append((data_gal[select][ind]['ra']<meanra)&(data_gal[select][ind]['dec']<meandec))
            patchgal.append((data_gal[select][ind]['ra']>meanra)&(data_gal[select][ind]['dec']<meandec))
            for pat in range(4):
                patchstarbool =  patchstars[pat]
                data_starsaux = data_stars[patchstarbool]
                patchgalbool = patchgal[pat]
                do_tau_stats( data_gal[select][ind][patchgalbool], Rs,
                              data_starsaux, bands, tilings, outpath,
                              max_sep=300, sep_units='arcmin', name= today + 'mod_bin_' + str(bin_c+ 1) + '_' + str(bin_c + 1) + '_patch_' + str(pat + 1),
                              bandcombo=args.bandcombo, mod=args.mod,
                              shapenoise=args.sn)

                
    galkeys = ['ra','dec','e_1','e_2','R11','R22']
    data_gal =  read_h5(args.metacal_cat, 'catalog/metacal/unsheared',  galkeys )
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
    #added by Lucas: modified to to select_2p and 2m
    Rs = [R11s, R22s]
    data_galaxies =  data_galaxies[select]
    print("Total objects after masking",  len(data_galaxies))
    print("R11s=",R11s)
    print("R22s=",R22s)

    patchstars = [];patchgal = []
    patchstars.append((data_stars['ra']>meanra)&(data_stars['dec']>meandec))
    patchstars.append((data_stars['ra']<meanra)&(data_stars['dec']>meandec))
    patchstars.append((data_stars['ra']<meanra)&(data_stars['dec']<meandec))
    patchstars.append((data_stars['ra']>meanra)&(data_stars['dec']<meandec))
    patchgal.append((data_gal[select]['ra']>meanra)&(data_galaxies[select]['dec']>meandec))
    patchgal.append((data_gal[select]['ra']<meanra)&(data_galaxies[select]['dec']>meandec))
    patchgal.append((data_gal[select]['ra']<meanra)&(data_galaxies[select]['dec']<meandec))
    patchgal.append((data_gal[select]['ra']>meanra)&(data_galaxies[select]['dec']<meandec))
    for pat in range(4):
        patchstarbool =  patchstars[pat]
        data_starsaux = data_stars[patchstarbool]
        patchgalbool = patchgal[pat]
        do_tau_stats( data_gal[patchgalbool], Rs,
                      data_starsaux, bands, tilings, outpath,
                      max_sep=300, sep_units='arcmin', name= today +
                      'mod_bin_' + str(bin_c+ 1) + '_patch_' + str(pat + 1),
                      bandcombo=args.bandcombo, mod=args.mod,
                      shapenoise=args.sn)

if __name__ == "__main__":
    main()
