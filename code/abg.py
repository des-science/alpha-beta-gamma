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
    from plot_stats import pretty_rho1, pretty_rho2, pretty_rho0,  pretty_tau
    from chi2 import CHI2,  minimizeCHI2,  plotCHI2, plotCHI2shifted
    from maxlikelihood import MCMC, OneParMaxLike
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
    plt.savefig(outpath +'/rho1_all_rsrs.pdf')
    plt.clf()
    pretty_rho2(meanr, rho2p, sig_rho2, sqrtn, rho5p, sig_rho5)
    plt.savefig(outpath +'/rho2_all_rsrs.pdf')
    plt.clf()
    pretty_rho0(meanr, rho0p, sig_rho0, sqrtn)
    plt.savefig(outpath +'/rho0_all_rsrs.pdf')

    #Reading and plotting reserved stars galaxies correlations
    meanr2, tau0p, tau2p, tau5p, sig_tau0, sig_tau2, sig_tau5 =  read_taus(args.rsgcorr)
    plt.clf()
    pretty_tau(meanr2, tau0p, sig_tau0, sqrtn, r'$\tau_{0}(\theta)$')
    plt.savefig(outpath +'/tau0_all_rsgal.pdf')
    plt.clf()
    pretty_tau(meanr2, tau2p, sig_tau2, sqrtn, r'$\tau_{2}(\theta)$')
    plt.savefig(outpath +'/tau2_all_rsgal.pdf')
    plt.clf()
    pretty_tau(meanr2, tau5p, sig_tau5, sqrtn, r'$\tau_{5}(\theta)$')
    plt.savefig(outpath +'/tau5_all_rsgal.pdf')
    
    
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

    #############################################################################
    #Residual Test
    '''
    gflag, bflag = False, False
    i_guess = [0] #fiducial values
    alpha0, chisq0 =  minimizeCHI2(data, i_guess,  eq=0, gflag=gflag, bflag=bflag)
    alpha1, chisq1 =  minimizeCHI2(data, i_guess,  eq=1, gflag=gflag, bflag=bflag)
    alpha2, chisq2 =  minimizeCHI2(data, i_guess,  eq=2, gflag=gflag, bflag=bflag)
    print(alpha0, alpha1, alpha2)
    
    res0 = (tau0p - alpha0*rho0p)**2 
    res1 = (tau2p - alpha1*rho2p)**2 
    res2 = (tau5p - alpha2*rho5p)**2
    plt.clf()
    plt.plot(meanr, sig_tau0**2, color='blue', label=r'$var(\tau_{0})$', marker='o')
    plt.plot(meanr, sig_tau2**2, color='red', label=r'$var(\tau_{2})$', marker='o')
    plt.plot(meanr, sig_tau5**2, color='green', label=r'$var(\tau_{5})$', marker='o')
    plt.legend(loc='lower left',  shadow=True,  fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    #plt.ylim( [10**-22,10 **-10] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel('Variances', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()
    print("Printing variances_bybin.pdf")
    plt.savefig(outpath +'/variances_bybin.pdf')
    plt.clf()
    plt.plot(meanr, res0, color='blue', label=r'$(\tau_{0}-\alpha_{0}\rho_{0})^2$', marker='o')
    plt.plot(meanr, res1, color='red', label=r'$(\tau_{2}-\alpha_{1}\rho_{2})^2$', marker='o')
    plt.plot(meanr, res2, color='green', label=r'$(\tau_{5}-\alpha_{2}\rho_{5})^2$', marker='o')
    plt.legend(loc='lower left',  shadow=True,  fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    plt.ylim( [10**-22,10 **-10] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel('Residuals', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()
    print("Printing Residual_chi2_bybin.pdf")
    plt.savefig(outpath +'/Residuals_chi2_bybin.pdf')
    plt.clf()
    plt.plot(meanr, res0/sig_tau0**2 , color='blue', label=r'$\chi_{0}^2$', marker='o')
    plt.plot(meanr, res1/sig_tau2**2, color='red', label=r'$\chi_{1}^2$', marker='o')
    plt.plot(meanr, res2/sig_tau5**2, color='green', label=r'$\chi_{2}^2$', marker='o')
    plt.legend(loc='upper left',  shadow=True,  fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    plt.ylim( [0.01,50000.] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\chi^{2}$', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()
    print("Printing chi2_bybin.pdf")
    plt.savefig(outpath +'/chi2_bybin.pdf')
    '''
    #############################################################################

    
    eq = 2
    nwalkers,  nsteps = 100,  1000
    dof = len(rhos[0])

    ## ALPHA
    gflag, bflag = False, False
    i_guess = [0] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    print("alpha" , fitted_params)
    print("reduced Chi2 :", chisq/dof )
    namemc = outpath+'/mcmc_alpha.pdf'
    namecont = outpath+'/contours_alpha.pdf'
    #MCMC(fitted_params,data,nwalkers,nsteps, namemc, namecont, eq=eq, gflag=gflag, bflag=bflag )
    '''
    x_arr= -0.05, 0.05, 100
    filename1 = outpath+'/chisq_only_alpha.pdf'
    filename2 = outpath+'/chisqshifted_only_alpha.pdf'
    plotCHI2(None, data,x_arr,filename1, eq=eq,gflag=gflag,bflag=bflag)
    plotCHI2shifted(None, data, x_arr, chisq, filename2, eq=eq, gflag=gflag, bflag=bflag)
    '''
    
    
    ## ALPHA-BETA
    gflag, bflag = False, True
    i_guess = [0,-1] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    print("alpha, beta" , fitted_params)
    print("reduced Chi2:", chisq/dof )
    namemc = outpath+'/mcmc_alpha-beta.pdf'
    namecont = outpath+'/contours_alpha-beta.pdf'
    #MCMC(fitted_params,data,nwalkers,nsteps, namemc, namecont, eq=eq, gflag=gflag, bflag=bflag )

    ## ALPHA-BETA-GAMMA
    gflag, bflag = True, True
    i_guess = [0,-1,- 1] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    print("alpha, beta, gamma:" , fitted_params)
    print("reduced Chi2:", chisq/dof )
    namemc = outpath+'/mcmc_alpha-beta-eta.pdf'
    namecont = outpath+'/contours_alpha-beta-eta.pdf'
    MCMC(fitted_params,data, nwalkers, nsteps, namemc, namecont,  eq=eq, gflag=gflag, bflag=bflag )

    
    #############################################################################
    #Decomposing by components
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
    #############################################################################
   

    
if __name__ == "__main__":
    main()
