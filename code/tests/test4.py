#Comparing rhostats
import os
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Alpha beta gamma test solving the system of equatiosn and plotting correlations')
    
    parser.add_argument('--rsgcorr',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_irz.json',
                        #default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/used/tau_all_galaxy-reserved_irz_newselection.json',
                        help='Json file with the reserved stars - galaxies correlations')
    parser.add_argument('--rsrscorr',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_mod_irz.json',
                        help='Json file with the reserved stars - reserved stars correlations')
    parser.add_argument('--outpath', default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/plots',
                        help='location of the output of the files')
    
    
    args = parser.parse_args()

    return args


def main():
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    from readjson import read_rhos, read_taus
    from plot_stats import pretty_rho

    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    #Reading a ploting reserved stars correlations
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 = read_rhos(args.rsrscorr)
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, rho2p + rho1p, sig_rho1, sqrtn, legend=r'$\rho_2(\theta)+\rho_1(\theta)$')
    print("Printing :", outpath +'/rho2+rho1_all_rsrs.pdf')
    plt.savefig(outpath +'/rho2+rho1_all_rsrs.pdf')


    usualrhos = "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_unmod_irz.json"
    modrhos =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_mod_irz.json"
    umeanr, urho0p, urho1p, urho2p, urho3p, urho4p, urho5p, usig_rho0, usig_rho1, usig_rho2, usig_rho3, usig_rho4, usig_rho5 = read_rhos(usualrhos)
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 = read_rhos(modrhos)
    plt.clf()
    pretty_rho(meanr, rho2p + rho1p, sig_rho1, sqrtn, legend=r'$\rho_2(\theta)+\rho_1(\theta)$')
    print("Printing :", outpath +'/rhoswithandwithoutmod.pdf')
    plt.savefig(outpath +'/rhoswithandwithoutmod.pdf')

if __name__ == "__main__":
    main()
   
