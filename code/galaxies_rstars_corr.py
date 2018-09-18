import os

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation galaxies and reserved stars')
    
    parser.add_argument('--metacal_cat',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3_master/Y3_mastercat_v2_6_20_18_subsampled.h5',
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--piff_cat',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3a1-v29',
                        help='Full Path to the Only stars Piff catalog')
    parser.add_argument('--exps_file',
                        default='/home/dfa/sobreira/alsina/DESWL/psf/ally3.grizY',
                        #default='/home/dfa/sobreira/alsina/DESWL/psf/testexp',
                        help='list of exposures (in lieu of separate exps)')
    parser.add_argument('--bands', default='grizY', type=str,
                         help='Limit to the given bands')
    parser.add_argument('--bandcombo', default=False,
                        action='store_const', const=True,
                        help='run rho2 for all combination of bands, if false run particular combination defined in band')
    parser.add_argument('--use_reserved', default=True,
                        action='store_const', const=True,
                        help='just use the objects with the RESERVED flag')
    parser.add_argument('--frac', default=1., type=float,
                        help='Choose a random fraction of the input stars')
    parser.add_argument('--outpath', default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma',
                        help='location of the output of the files')
    
    
    args = parser.parse_args()

    return args

        
def main():
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    import numpy as np
    import h5py as h
    from read_psf_cats import read_data
    from read_psf_cats import toList
    from run_rho import do_canonical_stats
    from run_rho import do_cross_stats
    
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
            'piff_e1', 'piff_e2', 'piff_T']
 
    exps = toList(args.exps_file)
    data_stars, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved)
    

    #Reading metacal catalog
    f = h.File(args.metacal_cat, 'r')
    unsheared =  f['catalog/metacal/unsheared']
    print('Unsheared Metacal catalog fields',unsheared.keys())
    data_galaxies =  {}
    data_galaxies['ra'] = np.array(unsheared['ra'])
    data_galaxies['dec'] = np.array(unsheared['dec'])
    data_galaxies['obs_e1'] = np.array(unsheared['e_1'])
    data_galaxies['obs_e2'] = np.array(unsheared['e_2'])
    #data_galaxies['band'] = np.array(unsheared['nimage_tot_i'])

    do_cross_stats(data_stars, data_galaxies, bands, tilings, outpath,
                      name='all_galaxy-reserved', bandcombo=args.bandcombo)
    

if __name__ == "__main__":
    main()
