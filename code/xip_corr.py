import os
today = '26-03-19_'

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Shear correlation function with metacal. ')
    
    parser.add_argument('--metacal_cat',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3_master/Y3_mastercat_v2_6_20_18_subsampled.h5',
                        #default='/home2/dfa/sobreira/alsina/catalogs/y3_master/cats_des_y3/Y3_mastercat_v2_6_20_18.h5', 
                        help='Full Path to the Metacalibration catalog')
    parser.add_argument('--bands', default='riz', type=str,
                         help='Limit to the given bands')
    parser.add_argument('--bandcombo', default=False,
                        action='store_const', const=True,
                        help='run rho2 for all combination of bands, if false run particular combination defined in bands')
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

def bin_pairs(nbins):
    import itertools
    a=[i for i in range(nbins)]
    bin_pairs=[]
    #for p in itertools.combinations_with_replacement(a,2):
    for p in itertools.product(a, repeat=2):
        bin_pairs.append(p)
    return bin_pairs 

def main():
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    #sys.path.insert(0, '/global/cscratch1/sd/alsina/alpha-beta-gamma/code/src')
    
    import numpy as np
    from read_psf_cats import read_h5
    from run_rho import  do_xi_stats,  do_xi_stats_tomo
    import h5py as h
    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
        
    if(args.tomo):
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
        
        ind2 = np.where( zbin_array==1 )[0]
        ind3 = np.where( zbin_array==2 )[0]
        ind4 = np.where( zbin_array==3 )[0]

        data_gal_list = []; Rs_list = []
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
            data_gal=  data_gal[select][ind]
            data_gal_list.append(data_gal)
            Rs_list.append(Rs)
        for i, j in bin_pairs(4):
            do_xi_stats_tomo(data_gal_list, Rs_list, i, j, outpath,
                             name= today +'xi_mod_'+ str(i + 1) +'_'
                             +str(j + 1) , max_sep=300,
                             sep_units='arcmin',
                             bandcombo=args.bandcombo, mod=args.mod,
                             shapenoise=args.sn)
    
    galkeys = ['ra','dec','e_1','e_2','R11','R22']
    data_galaxies =  read_h5(args.metacal_cat, 'catalog/metacal/unsheared',  galkeys )
    print("Total objects in catalog:", len(data_galaxies))
    dgamma = 2*0.01
    f = h.File(args.metacal_cat, 'r')
    index =  f['index']
    select = np.array(index['select'])
    select_1p = np.array(index['select_1p'])
    select_1m = np.array(index['select_1m'])
    select_2p = np.array(index['select_2p']) 
    select_2m = np.array(index['select_2m'])
    R11s = (data_galaxies['e_1'][select_1p].mean() - data_galaxies['e_1'][select_1m].mean() )/dgamma
    R22s = (data_galaxies['e_2'][select_2p].mean() - data_galaxies['e_2'][select_2m].mean() )/dgamma 
    Rs = [R11s, R22s]
    data_galaxies =  data_galaxies[select]
    print("Total objects after masking",  len(data_galaxies))
    print("R11s=",R11s)
    print("R22s=",R22s)
    
    do_xi_stats(data_galaxies, Rs, outpath, name= today + 'xi_mod',
                max_sep=300, sep_units='arcmin',
                bandcombo=args.bandcombo, mod=args.mod,
                shapenoise=args.sn)


if __name__ == "__main__":
    main()
