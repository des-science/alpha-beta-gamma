#What is happening with tau0 at large scales
import os

def main():
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    from readjson import read_rhos, read_taus
    from plot_stats import pretty_rho
    import numpy as np

    outpath = os.path.expanduser('/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/plots')
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
      
    tausp1 =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_mod_patch1_irz.json"
    tausp2 =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_mod_patch2_irz.json"
    tausp3 =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_mod_patch3_irz.json"
    tausp4 =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_mod_patch4_irz.json"
    meanr1, tau0p1, tau2p1, tau5p1, sig_tau01, sig_tau21, sig_tau51  = read_taus(tausp1)
    meanr2, tau0p2, tau2p2, tau5p2, sig_tau02, sig_tau22, sig_tau52  = read_taus(tausp2)
    meanr3, tau0p3, tau2p3, tau5p3, sig_tau03, sig_tau23, sig_tau53  = read_taus(tausp3)
    meanr4, tau0p4, tau2p4, tau5p4, sig_tau04, sig_tau24, sig_tau54  = read_taus(tausp4)
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr1, tau0p1, sig_tau01, sqrtn, legend='P++', lfontsize=10,  color='blue', ylabel=r'$\tau_{0}$', ylim=False)
    pretty_rho(meanr2, tau0p2, sig_tau02, sqrtn, legend='P-+', lfontsize=10,  color='red', ylabel=r'$\tau_{0}$', ylim=False)
    pretty_rho(meanr3, tau0p3, sig_tau03, sqrtn, legend='P--', lfontsize=10,  color='green', ylabel=r'$\tau_{0}$', ylim=False)
    pretty_rho(meanr4, tau0p4, sig_tau04, sqrtn, legend='P+-', lfontsize=10,  color='black', ylabel=r'$\tau_{0}$', ylim=False)
    print("Printing :", outpath +'/taus0_quadrants.png')
    plt.savefig(outpath +'/taus0_quadrants.png')

    plt.clf()
    pretty_rho(meanr1, tau2p1, sig_tau21, sqrtn, legend='P++', lfontsize=10,  color='blue', ylabel=r'$\tau_{2}$', ylim=False)
    pretty_rho(meanr2, tau2p2, sig_tau22, sqrtn, legend='P-+', lfontsize=10,  color='red', ylabel=r'$\tau_{2}$', ylim=False)
    pretty_rho(meanr3, tau2p3, sig_tau23, sqrtn, legend='P--', lfontsize=10,  color='green', ylabel=r'$\tau_{2}$', ylim=False)
    pretty_rho(meanr4, tau2p4, sig_tau24, sqrtn, legend='P+-', lfontsize=10,  color='black', ylabel=r'$\tau_{2}$', ylim=False)
    print("Printing :", outpath +'/taus2_quadrants.png')
    plt.savefig(outpath +'/taus2_quadrants.png')

    plt.clf()
    pretty_rho(meanr1, tau5p1, sig_tau51, sqrtn, legend='P++', lfontsize=10,  color='blue', ylabel=r'$\tau_{5}$', ylim=False)
    pretty_rho(meanr2, tau5p2, sig_tau52, sqrtn, legend='P-+', lfontsize=10,  color='red', ylabel=r'$\tau_{5}$', ylim=False)
    pretty_rho(meanr3, tau5p3, sig_tau53, sqrtn, legend='P--', lfontsize=10,  color='green', ylabel=r'$\tau_{5}$', ylim=False)
    pretty_rho(meanr4, tau5p4, sig_tau54, sqrtn, legend='P+-', lfontsize=10,  color='black', ylabel=r'$\tau_{5}$', ylim=False)
    print("Printing :", outpath +'/taus5_quadrants.png')
    plt.savefig(outpath +'/taus5_quadrants.png')


if __name__ == "__main__":
    main()
