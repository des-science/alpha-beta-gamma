import os

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Alpha beta gamma test solving the system of equatiosn and plotting correlations')
    
    parser.add_argument('--rsgcorr',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/sigma_all_galaxy-reserved_irz.json',
                        help='Json file with the reserved stars - galaxies correlations')
    parser.add_argument('--rsrscorr',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_irz.json',
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
    from readjson import read_rhos, read_sigmas
    from plot_stats import pretty_rho1, pretty_rho2, pretty_rho0,  pretty_sigma
    import numpy as np

    
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
    pretty_rho1(meanr, rho1p, sig_rho1, sqrtn, rho3p, sig_rho3, rho4p, sig_rho4)
    plt.savefig(outpath +'/rho1_all_rsrs.png')
    plt.clf()
    pretty_rho2(meanr, rho2p, sig_rho2, sqrtn, rho5p, sig_rho5)
    plt.savefig(outpath +'/rho2_all_rsrs.png')
    plt.clf()
    pretty_rho0(meanr, rho0p, sig_rho0, sqrtn)
    plt.savefig(outpath +'/rho0_all_rsrs.png')

    #Reading and plotting reserved stars galaxies correlations
    meanr2, sigma0p, sigma2p, sigma5p, sig_sigma0, sig_sigma2, sig_sigma5 =  read_sigmas(args.rsgcorr)
    plt.clf()
    pretty_sigma(meanr2, sigma0p, sig_sigma0, sqrtn, r'$\sigma_{0}(\theta)$')
    plt.savefig(outpath +'/sigma0_all_rsgal.png')
    plt.clf()
    pretty_sigma(meanr2, sigma2p, sig_sigma2, sqrtn, r'$\sigma_{2}(\theta)$')
    plt.savefig(outpath +'/sigma2_all_rsgal.png')
    plt.clf()
    pretty_sigma(meanr2, sigma5p, sig_sigma5, sqrtn, r'$\sigma_{5}(\theta)$')
    plt.savefig(outpath +'/sigma5_all_rsgal.png')

    #Solving equation system (AX=B) and ploting solution for each angular bin
    alpha = []
    beta = []
    gamma =  []
    for i in range(0, len(meanr)):
        A = np.array([[rho0p[i], rho2p[i], rho5p[i]], [rho2p[i], rho1p[i], rho4p[i]], [rho5p[i], rho4p[i], rho3p[i]]])
        B = np.array([sigma0p[i], sigma2p[i],  sigma5p[i]])
        sol = np.linalg.solve(A, B)
        alpha.append(sol[0])
        beta.append(sol[1])
        gamma.append(sol[2])
    #print(alpha, beta, gamma)
    plt.clf()
    plt.plot(meanr, alpha,  color='blue', label=r'$\alpha$', marker='o')
    plt.plot(meanr, beta,  color='red', label=r'$\beta$', marker='o')
    plt.plot(meanr, gamma,  color='green', label=r'$\gamma$', marker='o')
    plt.legend(loc='upper right', fontsize=24)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\alpha$,$\beta$,$\gamma$', fontsize=24)
    plt.xscale('log')
    #plt.yscale('log', nonposy='clip')
    plt.tight_layout()
    plt.savefig(outpath +'/alpha_beta_gamma.png')

if __name__ == "__main__":
    main()
