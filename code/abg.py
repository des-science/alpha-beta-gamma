import os

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Alpha beta gamma test solving the fitting problem of system of equatiosn, plotting correlations and final correlation function with bias')
    
    parser.add_argument('--rsgcorr',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_mod_irz.json',
                        help='Json file with the reserved stars - galaxies correlations')
    parser.add_argument('--rsrscorr',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_mod_epiff_magcut_sn_irz.json',
                        help='Json file with the reserved stars - reserved stars correlations')
    parser.add_argument('--xipobs',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/xi_xi_mod_riz.json',
                        help='Json file with the xip_obs correlation function')
    parser.add_argument('--outpath', default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/plots',
                        help='location of the output of the files')
    parser.add_argument('--maxscale', default=None, type=float, 
                        help='Limit the analysis to certain maximum scale, units are determined by .json file with the correlations')
    
    
    args = parser.parse_args()

    return args

def plotallrhos(args, outpath):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    from readjson import read_rhos
    from plot_stats import pretty_rho1, pretty_rho2, pretty_rho0,  pretty_tau
    
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 = read_rhos(args.rsrscorr, args.maxscale)
    sqrtn = 1
    plt.clf()
    pretty_rho1(meanr, rho1p, sig_rho1, sqrtn, rho3p, sig_rho3, rho4p, sig_rho4)
    print("Printing file: ", outpath +'/rho1_all_rsrs.png')
    plt.savefig(outpath +'/rho1_all_rsrs.png')
    plt.clf()
    pretty_rho2(meanr, rho2p, sig_rho2, sqrtn, rho5p, sig_rho5)
    print("Printing file: ", outpath +'/rho2_all_rsrs.png')
    plt.savefig(outpath +'/rho2_all_rsrs.png')
    plt.clf()
    pretty_rho0(meanr, rho0p, sig_rho0, sqrtn)
    print("Printing file: ", outpath +'/rho0_all_rsrs.png')
    plt.savefig(outpath +'/rho0_all_rsrs.png')
    
def plotalltaus(args, outpath):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    from readjson import read_taus
    from plot_stats import pretty_tau
    
    meanr2, tau0p, tau2p, tau5p, sig_tau0, sig_tau2, sig_tau5 =  read_taus(args.rsgcorr, args.maxscale)
    sqrtn = 1
    plt.clf()
    pretty_tau(meanr2, tau0p, sig_tau0, sqrtn, r'$\tau_{0}(\theta)$')
    print("Printing file: ", outpath +'/tau0_all_rsgal.png')
    plt.savefig(outpath +'/tau0_all_rsgal.png')
    plt.clf()
    pretty_tau(meanr2, tau2p, sig_tau2, sqrtn, r'$\tau_{2}(\theta)$')
    print("Printing file: ", outpath +'/tau2_all_rsgal.png')
    plt.savefig(outpath +'/tau2_all_rsgal.png')
    plt.clf()
    pretty_tau(meanr2, tau5p, sig_tau5, sqrtn, r'$\tau_{5}(\theta)$')
    print("Printing file: ", outpath +'/tau5_all_rsgal.png')
    plt.savefig(outpath +'/tau5_all_rsgal.png')

    
def plotxipandbias(samples, rhosfilename, tausfilename, xifilename, outpath):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    from readjson import read_rhos, read_taus,  read_xi
    from plot_stats import pretty_rho
    from maxlikelihood import bestparameters
    import numpy as np

    a, b, n = bestparameters(samples)
    print(a, b, n)
    par_matcov = np.cov(samples)
    #print(par_matcov)
    vara, varb, varn = par_matcov[0, 0], par_matcov[1, 1], par_matcov[2, 2]
    #print(np.sqrt(vara), np.sqrt(varb), np.sqrt(varn))
    covab, covan, covbn = par_matcov[0, 1], par_matcov[0, 2], par_matcov[1, 2]
    
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p,\
            sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4,\
            sig_rho5 = read_rhos(rhosfilename)
    meanr, tau0p, tau2p, tau5p, sig_tau0, sig_tau2,sig_tau5\
            = read_taus(tausfilename)
    meanr, xip_obs,  sig_xip = read_xi(xifilename)

    #Ploting each term of the bias
    
    sqrtn = 1
    #supposing that a,b and n are idependent of rhos(scale independent)
    var0 = ((2*a*rho0p)**2)*vara +  (a**2)*(sig_rho0**2)
    var1 = ((2*b*rho1p)**2)*varb +  (b**2)*(sig_rho1**2)
    var2 = ((2*n*rho3p)**2)*varn +  (n**2)*(sig_rho3**2)
    varab = ((a*b)**2)*( (vara/((a)**2)) + (varb/((b)**2)) + 2*covab/(a*b) ) 
    #print(varab)
    var3 = 4*((a*b*rho2p)**2)*(varab/((a*b)**2) + (sig_rho2/rho2p)**2)
    varbn = ((n*b)**2)*( (varn/((n)**2)) + (varb/((b)**2)) + 2*covbn/(b*n) )
    #print(varbn)
    var4 = 4*((n*b*rho4p)**2)*(varbn/((b*n)**2) + (sig_rho4/rho4p)**2)
    varan = ((n*a)**2)*( (varn/((n)**2)) + (vara/((a)**2)) + 2*covan/(a*n) ) 
    #print(varan)
    var5 = 4*((n*a*rho5p)**2)*(varan/((a*n)**2) + (sig_rho5/rho5p)**2) 
    plt.clf()
    
    pretty_rho(meanr, (a**2)*rho0p, np.sqrt(var0), sqrtn, legend=r'$\alpha^{2} \rho_{0}$',lfontsize=10,  color='red', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, (b**2)*rho1p, np.sqrt(var1), sqrtn, legend=r'$\beta^{2}\rho_{1}$',lfontsize=10,  color='green', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, (n**2)*rho3p, np.sqrt(var2), sqrtn, legend=r'$\eta^{2}\rho_{3}$', lfontsize=10, color='black', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, (2*a*b)*rho2p, np.sqrt(var3), sqrtn, legend=r'$2\alpha\beta \rho_{2}$',lfontsize=10,  color='yellow', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, (2*b*n)*rho4p, np.sqrt(var4), sqrtn, legend=r'$2\beta\eta\rho_{4}$',lfontsize=10,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, (2*n*a)*rho5p, np.sqrt(var5), sqrtn, legend=r'$2\eta\alpha\rho_{5}$', lfontsize=10, color='gray', ylabel='Correlations', ylim=False)

    fname = '/xibias_parts.pdf' 
    print(outpath +fname)
    plt.savefig(outpath +fname)
    
    dxip = (a**2)*rho0p + (b**2)*rho1p + (n**2)*rho3p + (2*a*b)*rho2p + (2*b*n)*rho4p + (2*n*a)*rho5p
    f1 = 2*(a*rho0p + b*rho2p + n*rho5p)     
    f2 = 2*(b*rho1p + a*rho2p + n*rho4p)
    f3 = 2*(n*rho3p + b*rho4p + a*rho5p)
    f4 = a**2 ; f5 = b**2; f6 = 2*a*b
    f7 = n**2 ; f8 = 2*b*n; f9 = 2*n*a 
    var_dxip = (f1**2)*vara + (f2**2)*varb + (f3**2)*varn + \
    (f4**2)*(sig_rho0**2) + (f5**2)*(sig_rho1**2) + \
    (f6**2)*(sig_rho2**2) + (f7**2)*(sig_rho3**2) + \
    (f8**2)*(sig_rho4**2) + (f9**2)*(sig_rho5**2) + 2*(f1*f2*covab + f1*f3*covan + f2*f3*covbn)         
    var_dxip2 = var0 + var1 + var2 + var3 + var4 + var5
    print(var_dxip)
    print(var_dxip2)
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, xip_obs, sig_xip, sqrtn, legend=r"$\xi_{+}^{obs}$",lfontsize=10,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, dxip, np.sqrt(var_dxip), sqrtn, legend=r"$\delta \xi_{+}$",lfontsize=10,  color='red', ylabel='Correlations', ylim=False)
    print(sig_xip)
    print(np.sqrt(var_dxip))
    fname = '/xiobs_vs_xibias.pdf' 
    print(outpath +fname)
    plt.savefig(outpath +fname)

    
def main():
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    from readjson import read_rhos, read_taus
    from chi2 import CHI2,  minimizeCHI2,  plotCHI2, plotCHI2shifted
    from maxlikelihood import MCMC, percentiles, OneParMaxLike
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
    plotallrhos(args, outpath=outpath)

    #Reading and plotting reserved stars galaxies correlations
    plotalltaus(args, outpath=outpath)
   
    
    
    #Finding best alpha beta gamma
    
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 = read_rhos(args.rsrscorr, args.maxscale)
    rhos = [rho0p, rho1p, rho2p, rho3p, rho4p, rho5p]
    sigrhos = [sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5]
    meanr2, tau0p, tau2p, tau5p, sig_tau0, sig_tau2, sig_tau5 =  read_taus(args.rsgcorr, args.maxscale)
    taus = [tau0p, tau2p, tau5p]
    sigtaus = [sig_tau0, sig_tau2, sig_tau5]
    data = {}
    data['rhos'] = rhos
    data['sigrhos'] = sigrhos
    data['taus'] = taus
    data['sigtaus'] = sigtaus
    
    
    nwalkers,  nsteps = 100,  1000
    dof = len(rhos[0])
    moderr = False
    nsig = 1

    #for eq in [0, 1, 2, 4]:
    for eq in [None]:
        '''
        ## ALPHA
        gflag, bflag = False, False
        i_guess = [0] #fiducial values
        fitted_params, chisq = minimizeCHI2(data, i_guess, eq=eq,
                                            gflag=gflag, bflag=bflag,
                                            moderr=moderr)
        print("alpha" , fitted_params)
        print("reduced Chi2 :", chisq/dof )
        namemc = None
        #namemc = outpath+'/mcmc_alpha_eq' + str(eq) + '_.pdf'
        namecont = outpath+'/contours_alpha_eq' + str(eq) + '_.pdf'
        samples = MCMC(fitted_params,data,nwalkers,nsteps, namemc,
                       namecont, eq=eq, gflag=gflag, bflag=bflag,
                       moderr=moderr )
        mcmcpars = percentiles(samples, nsig=nsig) 
        print("mcmc_alpha",  mcmcpars)    
    
        ## ALPHA-BETA
        gflag, bflag = False, True
        i_guess = [0,-1] #fiducial values
        fitted_params, chisq = minimizeCHI2(data, i_guess, eq=eq,
                                            gflag=gflag, bflag=bflag,
                                            moderr=moderr)
        print("alpha, beta" , fitted_params)
        print("reduced Chi2:", chisq/dof )
        namemc =  None
        #namemc = outpath+'/mcmc_alpha-beta_eq' + str(eq) + '_.pdf'
        namecont = outpath+'/contours_alpha-beta_eq' + str(eq) + '_.pdf'
        samples = MCMC(fitted_params,data,nwalkers,nsteps, namemc,
                       namecont, eq=eq, gflag=gflag, bflag=bflag,
                       moderr=moderr )
        mcmcpars = percentiles(samples, nsig=nsig) 
        print("mcmc_alpha-beta",  mcmcpars)    
        '''
        
        ## ALPHA-BETA-ETA
        gflag, bflag = True, True
        i_guess = [0,-1,- 1] #fiducial values
        fitted_params, chisq = minimizeCHI2(data, i_guess, eq=eq,
                                            gflag=gflag, bflag=bflag,
                                            moderr=moderr)
        print("alpha, beta, gamma:" , fitted_params)
        print("reduced Chi2:", chisq/dof )
        namemc = outpath+'/mcmc_alpha-beta-eta_eq' + str(eq) + '_.png'
        namecont = outpath+'/contours_alpha-beta-eta_eq' + str(eq) + '_.png'
        samples = MCMC(fitted_params,data, nwalkers, nsteps, namemc,
                       namecont, eq=eq, gflag=gflag, bflag=bflag,
                       moderr=moderr,  plot=False )
        #print(samples)
        mcmcpars = percentiles(samples, nsig=nsig) 
        print("mcmc_alpha-beta-eta",  mcmcpars)
        
    plotxipandbias(samples, rhosfilename=args.rsrscorr , tausfilename=args.rsgcorr, xifilename=args.xipobs, outpath=outpath)
        
        
   
if __name__ == "__main__":
    main()
