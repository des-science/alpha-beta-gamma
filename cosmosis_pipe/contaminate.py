import os

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Sofware to plot all quantities after running WL-pipeline')
    
    parser.add_argument('--original',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/cosmosis_pipe/2pt_sim_1110_baseline_Y3cov.fits',
                        help='File containing xip to be modified')
    parser.add_argument('--contaminant',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/cosmosis_pipe/abg_dxip_tomo.fits',
                        help='Path for the outputs of this code')
    parser.add_argument('--filename',
                        default='2pt_sim_1110_baseline_Y3cov_contaminated.fits',
                        help='Path for the outputs of this code')  
    
    args = parser.parse_args()
    return args


  


def main():
    import fitsio
    from fitsio import FITS,FITSHDR
    args = parse_args()
    
    fitori = FITS(args.original)
    xiptable = fitori['xip']
    covmatrix = fitori['COVMAT']
    print(xiptable['BIN1'].read())

    dxiptable =  fitsio.read(args.contaminant,  ext=1)
    print(fitcontable)

    nbins=4
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin1 = (xiptable['BIN1']==i)&(xiptable['BIN2']==j)
        bin2 = (dxiptable['BIN1']==i)&(dxiptable['BIN2']==j)
        xiptable['VALUE'][bin1] -= dxiptable['VALUE'][bin2]
if __name__ == "__main__":
    main()
