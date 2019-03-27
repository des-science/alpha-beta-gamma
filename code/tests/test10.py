#Finding priors for beta and eta
import os
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Finding priors for beta and eta')

    parser.add_argument('--metacal_cat',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3_master/cats_des_y3/Y3_mastercat_v2_6_20_18.h5', 
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
    
    
    args = parser.parse_args()

    return args

def main():
    
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    import numpy as np
    from read_psf_cats import read_data, toList, read_h5
    import h5py as h
    args = parse_args()

    #Reading Metacal
    galkeys = ['ra', 'dec', 'T', 'psf_T']
    data_galaxies =  read_h5(args.metacal_cat, 'catalog/metacal/unsheared',  galkeys )
    print("Total objects in catalog:", len(data_galaxies))
    
    f = h.File(args.metacal_cat, 'r')
    index =  f['index']
    select = np.array(index['select'])
    data_galaxies =  data_galaxies[select]
    print("Total objects after masking",  len(data_galaxies))
    
    Tpsfgal = data_galaxies['psf_T']
    Tgal = data_galaxies['T']
    prior_gal = np.mean(Tpsfgal/Tgal)
    priorgal_max = np.max(Tpsfgal/Tgal)
    priorgal_min = np.min(Tpsfgal/Tgal)
    print('prior_gal=',  prior_gal)
    print('prior_gal_min=',  priorgal_min)
    print('prior_gal_max=',  priorgal_max)

    
    
    #Reading Mike stars catalog
    keys = ['ra', 'dec', 'obs_T', 'piff_T', 'mag']
 
    exps = toList(args.exps_file)
    data_stars, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved)
    print("Objects",  len(data_stars))
    data_stars = data_stars[data_stars['mag']<20]
    print("Objects with magnitude <20",  len(data_stars))
    Tpsf = data_stars['piff_T']
    Tstars = data_stars['obs_T']
    prior = np.mean(Tpsf/Tstars)
    prior_max = np.max(Tpsf/Tstars)
    prior_min = np.min(Tpsf/Tstars)
    print('prior_stars=',  prior)
    print('prior_stars_min=',  prior_min)
    print('prior_stars_max=',  prior_max)


    
    

if __name__ == "__main__":
    main()
    
