#plotting each term in the equation of correlations, to see if there are reason to cancel some part of the model.
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='')
    
    parser.add_argument('--rsgcorr',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_irz.json',
                        help='Json file with the reserved stars - galaxies correlations')
    parser.add_argument('--rsrscorr',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_mod_irz.json',
                        help='Json file with the reserved stars - reserved stars correlations')
    parser.add_argument('--outpath', default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/plots',
                        help='location of the output of the files')
    
    
    args = parser.parse_args()

    return args
def main():
    import os
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    from readjson import read_rhos, read_taus
    from chi2 import minimizeCHI2

    args = parse_args()
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

        
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 = read_rhos(args.rsrscorr)
    meanr2, tau0p, tau2p, tau5p, sig_tau0, sig_tau2, sig_tau5 =  read_taus(args.rsgcorr)
    
    #Finding best alpha beta gamma
    rhos = [rho0p, rho1p, rho2p, rho3p, rho4p, rho5p]
    sigrhos = [sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5]
    taus = [tau0p, tau2p, tau5p]
    sigtaus = [sig_tau0, sig_tau2, sig_tau5]
    data = {}
    data['rhos'] = rhos
    data['sigrhos'] = sigrhos
    data['taus'] = taus
    data['sigtaus'] = sigtaus

    gflag, bflag = True, True
    eq = None
    i_guess = [0,-1,- 1] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
  
    alpha, beta, eta = fitted_params
    plt.clf()
    plt.plot(meanr, tau0p, color='blue', label=r'$\tau_{0}$', marker='o')
    plt.plot(meanr, -tau0p, color='blue', ls=':', marker='o')
    plt.plot(meanr, alpha*rho0p, color='red', label=r'$\alpha \rho_{0}$', marker='o')
    plt.plot(meanr, -alpha*rho0p, color='red',ls=':', marker='o')
    plt.plot(meanr, beta*rho2p, color='green', label=r'$\beta\rho_{2}$', marker='o')
    plt.plot(meanr,-beta*rho2p, color='green',ls=':', marker='o')
    plt.plot(meanr, eta*rho5p, color='black', label=r'$\eta\rho_{5}$', marker='o')
    plt.plot(meanr,-eta*rho5p, color='black',ls=':', marker='o')
    plt.legend(loc='lower left',  shadow=True,  fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    #plt.ylim( [10**-22,10 **-10] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel('Correlation', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()
    print("Printing corr_pars1.pdf")
    plt.savefig(outpath +'/corr_pars1.pdf')
    
    plt.clf()
    plt.plot(meanr, tau2p, color='blue', label=r'$\tau_{2}$', marker='o')
    plt.plot(meanr, -tau2p, color='blue', ls=':', marker='o')
    plt.plot(meanr, alpha*rho2p, color='red', label=r'$\alpha \rho_{2}$', marker='o')
    plt.plot(meanr, -alpha*rho2p, color='red',ls=':', marker='o')
    plt.plot(meanr, beta*rho1p, color='green', label=r'$\beta\rho_{1}$', marker='o')
    plt.plot(meanr,-beta*rho1p, color='green',ls=':', marker='o')
    plt.plot(meanr, eta*rho4p, color='black', label=r'$\eta\rho_{4}$', marker='o')
    plt.plot(meanr,-eta*rho4p, color='black', ls=':',marker='o')
    plt.plot(meanr, alpha*rho2p + beta*rho1p + eta*rho4p, color='blue', label='rhs', marker='P')
    plt.plot(meanr,-alpha*rho2p -  beta*rho1p -  eta*rho4p, color='blue', ls=':', marker='P')
    plt.legend(loc='lower left',  shadow=True,  fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    #plt.ylim( [10**-22,10 **-10] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel('Correlation', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()
    print("Printing corr_pars2.pdf")
    plt.savefig(outpath +'/corr_pars2.pdf')

    plt.clf()
    plt.plot(meanr, tau5p, color='blue', label=r'$\tau_{5}$', marker='o')
    plt.plot(meanr, -tau5p, color='blue', ls=':', marker='o')
    plt.plot(meanr, alpha*rho5p, color='red', label=r'$\alpha \rho_{5}$', marker='o')
    plt.plot(meanr, -alpha*rho5p, color='red',ls=':', marker='o')
    plt.plot(meanr, beta*rho4p, color='green', label=r'$\beta\rho_{4}$', marker='o')
    plt.plot(meanr,-beta*rho4p, color='green',ls=':', marker='o')
    plt.plot(meanr, eta*rho3p, color='black', label=r'$\eta\rho_{3}$', marker='o')
    plt.plot(meanr,-eta*rho3p, color='black', ls=':', marker='o')
    plt.plot(meanr, alpha*rho5p + beta*rho4p + eta*rho3p, color='blue', label='rhs', marker='P')
    plt.plot(meanr,-alpha*rho5p -  beta*rho4p -  eta*rho3p, color='blue', ls=':', marker='P')
    plt.legend(loc='upper right',  shadow=True,  fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    #plt.ylim( [10**-22,10 **-10] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel('Correlation', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()
    print("Printing corr_pars3.pdf")
    plt.savefig(outpath +'/corr_pars3.pdf')
    
    
    
if __name__ == "__main__":
    main()
