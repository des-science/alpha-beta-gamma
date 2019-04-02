#Comparing bias with xiobs and parts
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Plot to compare xipobs with possible due to psf modelling ')
    
    parser.add_argument('--xipobs',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/xi_xi_mod_riz.json',
                        help='Json file with the xip_obs correlation function')
    parser.add_argument('--dxip',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/dxip-alpha-beta-eta.json',
                        help='Json file with the xip_obs correlation function')
    parser.add_argument('--srcpath',
                         default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src',
                         help='Path to src/ folder with plenty functions required')
    parser.add_argument('--outpath',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/',
                        help='location of the output of the files, by default plots is created')
    
    
    args = parser.parse_args()

    return args

def main():
    import sys
    args = parse_args()
    sys.path.insert(0, args.srcpath)
    import numpy as np
    from readjson import read_xi,  read_dxip
    from plot_stats import pretty_rho

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    
    meanr_xip, xip_obs,  sig_xip = read_xi(args.xipobs)
    meanr_dxip, dxip,  sig_dxip = read_dxip(args.dxip)
    plt.clf()
    pretty_rho(meanr_xip, xip_obs, sig_xip, legend=r"$\xi_{+}^{obs}$",lfontsize=10,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr_dxip, dxip, sig_dxip,  legend=r"$\delta \xi_{+}$",lfontsize=10,  color='red', ylabel='Correlations', ylim=False)
    fname = 'plots/xiobs_vs_xibias_abe2.png' 
    print('Printing:',  outpath +fname)
    plt.savefig(outpath +fname)

    meanr_dxip, dxip,  sig_dxip = read_dxip('/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/dxip-alpha-beta.json')
    plt.clf()
    pretty_rho(meanr_xip, xip_obs, sig_xip, legend=r"$\xi_{+}^{obs}$",lfontsize=10,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr_dxip, dxip, sig_dxip,  legend=r"$\delta \xi_{+}$",lfontsize=10,  color='red', ylabel='Correlations', ylim=False)
    fname = 'plots/xiobs_vs_xibias_ab2.png' 
    print('Printing:',  outpath +fname)
    plt.savefig(outpath +fname)

    meanr_dxip, dxip,  sig_dxip = read_dxip('/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/dxip-alpha.json')
    plt.clf()
    pretty_rho(meanr_xip, xip_obs, sig_xip, legend=r"$\xi_{+}^{obs}$",lfontsize=10,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr_dxip, dxip, sig_dxip,  legend=r"$\delta \xi_{+}$",lfontsize=10,  color='red', ylabel='Correlations', ylim=False)
    fname = 'plots/xiobs_vs_xibias_a2.png' 
    print('Printing:',  outpath +fname)
    plt.savefig(outpath +fname)
    
  


    
    

if __name__ == "__main__":
    main()
    
