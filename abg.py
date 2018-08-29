def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Read catalog program format h5')
    
    parser.add_argument('--psf_cat', default='/home2/dfa/sobreira/alsina/catalogs/y3_master/Y3_mastercat_v1_6_20_18_subsampled.h5',
                        help='Full Path to the catalog')
 

    args = parser.parse_args()

    return args


def main():
    from astropy.io import fits
    import fitsio
    import numpy as np
    import pandas
    import h5py as h
    
    args = parse_args()
    #psf_data=fitsio.read(args.psf_cat)
    f = h.File(args.psf_cat, 'r')
    print(f.keys())

if __name__ == "__main__":
    main()
