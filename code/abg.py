import os

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Alpha beta gamma test solving the system of equatiosn and plotting correlations')
    
    parser.add_argument('--rsgcorr',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/blabla.json',
                        help='Json file with the reserved stars - galaxies correlations')
    parser.add_argument('--rsrscorr',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/blabla.json',
                        help='Json file with the reserved stars - reserved stars correlations')
    parser.add_argument('--outpath', default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/plots',
                        help='location of the output of the files')
    
    
    args = parser.parse_args()

    return args

        
def main():
    import numpy as np

    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
        

    #STATISTIC USING ONLY RESERVED STARS
    keys = ['ra', 'dec','obs_e1', 'obs_e2', 'obs_T',
            'piff_e1', 'piff_e2', 'piff_T']
 
    exps = toList(args.exps_file)
    data, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved)
    do_canonical_stats(data, bands, tilings, outpath,
                       name='all_reserved', bandcombo=args.bandcombo)
    

    #CROSS STATISTIC GALAXY-STAR
    f = h.File(args.metacal_cat, 'r')
    unsheared =  f['catalog/metacal/unsheared']
    print('Unsheared Metacal catalog fields',unsheared.keys())
   
    

if __name__ == "__main__":
    main()
