import os

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Shear correlation function with metacal. ')
    
    parser.add_argument('--metacal_cat',
                        #default='/home2/dfa/sobreira/alsina/catalogs/y3_master/Y3_mastercat_v2_6_20_18_subsampled.h5',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3_master/cats_des_y3/Y3_mastercat_v2_6_20_18.h5', 
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--bands', default='riz', type=str,
                         help='Limit to the given bands')
    parser.add_argument('--bandcombo', default=False,
                        action='store_const', const=True,
                        help='run rho2 for all combination of bands, if false run particular combination defined in bands')
    parser.add_argument('--outpath', default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma',
                        help='location of the output of the files')
    
    
    args = parser.parse_args()

    return args

        
def main():
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    #sys.path.insert(0, '/global/cscratch1/sd/alsina/alpha-beta-gamma/code/src')
    
    import numpy as np
    from read_psf_cats import read_h5
    from run_rho import  do_xi_stats
    import h5py as h
    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
        

    
    galkeys = ['ra','dec','e_1','e_2','snr','size_ratio','flags','T','T_err','R11','R22']
    data_galaxies =  read_h5(args.metacal_cat, 'catalog/metacal/unsheared',  galkeys )
    print("Total objects in catalog:", len(data_galaxies))
    dgamma = 2*0.01
    f = h.File(args.metacal_cat, 'r')
    index =  f['index']
    select = np.array(index['select'])
    select_1p = np.array(index['select_1p'])
    select_1m = np.array(index['select_1m'])
    select_2p = np.array(index['select_2p']) #added by Lucas: each component gets sheared and this leads to different selection
    select_2m = np.array(index['select_2m']) #added by Lucas
    R11s = (data_galaxies['e_1'][select_1p].mean() - data_galaxies['e_1'][select_1m].mean() )/dgamma
    R22s = (data_galaxies['e_2'][select_2p].mean() - data_galaxies['e_2'][select_2m].mean() )/dgamma #added by Lucas: modified to to select_2p and 2m
    Rs = [R11s, R22s]
    data_galaxies =  data_galaxies[select]
    print("Total objects after masking",  len(data_galaxies))
    print("R11s=",R11s)
    print("R22s=",R22s)
    
    do_xi_stats(data_galaxies, Rs, outpath, name='xi_mod', bandcombo=args.bandcombo, mod=True, shapenoise=True)


if __name__ == "__main__":
    main()
