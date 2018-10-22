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

    tausepiff =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_max_sep30_irz.json"
    meanr, tau0p, tau2p, tau5p, sig_tau0, sig_tau2, sig_tau5  = read_taus(tausepiff)
    
    sqrtn = 1
    marker = 'o'
    plt.clf()
    plt.plot(meanr, tau0p, color='blue', label=r"$\tau_{0}$")
    plt.plot(meanr, -tau0p, color='blue', ls=':')
    plt.errorbar(meanr[tau0p>0], tau0p[tau0p>0], yerr=sig_tau0[tau0p>0]/sqrtn, color='blue', ls='', marker=marker)
    plt.errorbar(meanr[tau0p<0], -tau0p[tau0p<0], yerr=sig_tau0[tau0p<0]/sqrtn, color='blue', ls='', marker=marker)
    tau0p0_line = plt.errorbar(-meanr, tau0p, yerr=sig_tau0, color='blue', marker='o')

    plt.plot(meanr, tau2p, color='green', label=r"$\tau_{2}$")
    plt.plot(meanr, -tau2p, color='green', ls=':')
    plt.errorbar(meanr[tau2p>0], tau2p[tau2p>0], yerr=sig_tau2[tau2p>0]/sqrtn, color='green', ls='', marker=marker)
    plt.errorbar(meanr[tau2p<0], -tau2p[tau2p<0], yerr=sig_tau2[tau2p<0]/sqrtn, color='green', ls='', marker=marker)
    tau0p0_line = plt.errorbar(-meanr, tau0p, yerr=sig_tau0, color='green', marker='o')

    plt.plot(meanr, tau5p, color='red', label=r"$\tau_{5}$")
    plt.plot(meanr, -tau5p, color='red', ls=':')
    plt.errorbar(meanr[tau5p>0], tau5p[tau5p>0], yerr=sig_tau5[tau5p>0]/sqrtn, color='red', ls='', marker=marker)
    plt.errorbar(meanr[tau5p<0], -tau5p[tau5p<0], yerr=sig_tau5[tau5p<0]/sqrtn, color='red', ls='', marker=marker)
    tau0p0_line = plt.errorbar(-meanr, tau0p, yerr=sig_tau0, color='red', marker='o')

    plt.legend(loc='best', fontsize=24)
    plt.ylim( [1.e-11, 1.e-6] )
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,150.] )
    plt.xlabel(r'$\theta$ (degrees)', fontsize=24)
    plt.ylabel(r'$\tau$', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()
   
 
    print("Printing :", outpath +'/taus_30degrees.pdf')
    plt.savefig(outpath +'/taus_30degrees.pdf')

if __name__ == "__main__":
    main()
