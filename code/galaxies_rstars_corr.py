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
    from read_psf_cats import read_data, toList, read_h5
    from run_rho import do_canonical_stats,  do_cross_stats
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
            'piff_e1', 'piff_e2', 'piff_T']
 
    exps = toList(args.exps_file)
    data_stars, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved)
    
    
    galkeys = ['ra','dec','e_1','e_2','snr','size_ratio','flags','T','T_err','R11','R22']
    data_galaxies =  read_h5(args.metacal_cat, 'catalog/metacal/unsheared',  galkeys )

    dgamma = 2*0.01
    f = h.File(args.metacal_cat, 'r')
    index =  f['index']
    select = np.array(index['select'])
    select_1p = np.array(index['select_1p'])
    select_1m = np.array(index['select_1m'])
    R11s = (data_galaxies['e_1'][select_1p].mean() - data_galaxies['e_1'][select_1m].mean() )/dgamma
    R22s = (data_galaxies['e_2'][select_1p].mean() - data_galaxies['e_2'][select_1m].mean() )/dgamma
    Rs = [R11s, R22s]
    data_galaxies =  data_galaxies[select]
    
    '''
    print(len(data_galaxies))
    mask =  (data_galaxies['snr'] > 10)
    mask &= (data_galaxies['snr'] < 100)
    mask &= (data_galaxies['size_ratio']>0.5)
    mask &= (data_galaxies['flags']==0)
    mask &= ( (data_galaxies['T'] - 2 * data_galaxies['T_err']) > 0)
    data_galaxies =  data_galaxies[mask]
    print(len(data_galaxies))
    '''
    do_cross_stats(data_stars, data_galaxies, Rs, bands, tilings, outpath,
                   name='all_galaxy-reserved', bandcombo=args.bandcombo, mod=True)


if __name__ == "__main__":
    main()
