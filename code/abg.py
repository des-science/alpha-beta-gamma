import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Alpha beta gamma test solving the fitting problem of system ofequatiosn, plotting correlations and final correlation function withbias')
    
    parser.add_argument('--taus',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tomo_taus/tau_irz_01-04-19_all_galaxy-reserved_mod.json',
                        help='Json file with the reserved stars -galaxies correlations')
    parser.add_argument('--rhos',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_irz_26-03-19_all_reserved_mod_epiff_magcut_sn.json',
                        help='Json file with the reserved stars -reserved stars correlations')
    parser.add_argument('--maxscale', default=None, type=float, 
                        help='Limit the analysis to certain maximum scale, units are determined by .json file with the correlations')
    parser.add_argument('--uwmprior', default=False,
                        action='store_const', const=True, help='Run all tomographic correlations')
    parser.add_argument('--plots', default=False,
                        action='store_const', const=True, help='Plot correlations functions')
    parser.add_argument('--outpath',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/',
                        help='location of the output of the files, by default plots is created')
    parser.add_argument('--srcpath',
                         default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src',
                         help='Path to src/ folder with plenty functions required')
    parser.add_argument('--abe', default=False,
                        action='store_const', const=True, help='Run alpha, beta and eta test.')
    parser.add_argument('--ab', default=False,
                        action='store_const', const=True, help='Run alpha and beta test.')
    parser.add_argument('--a', default=False,
                        action='store_const', const=True, help='Run only alpha test.')
    
    
    
    args = parser.parse_args()

    return args
def writedxip(meanr, dxip, sig_dxip,  filename='dxip.json' ):
    import json
    stats = [meanr.tolist(), dxip.tolist(),  sig_dxip.tolist()]
    with open(filename,'w') as fp:
        json.dump([stats], fp)
    print('Done writting',  filename)
def writedxip_txt(meanr, dxip, sig_dxip, outpath,  namer='meanr.txt', namedxip='dxip.txt', namesdxip ='sig_dxip.txt' ):
    import numpy as np
    np.savetxt(outpath + namer, meanr,  header='angular separation in arcmin')
    np.savetxt(outpath + namedxip, dxip,  header='delta xip: bias of xipobs fue to psf modeling')
    np.savetxt(outpath + namesdxip, sig_dxip,  header='error = one standard deviation')

def getxipbias(samples, rhosfilename, plotname='terms_dxip.png',  plots=False):
    from readjson import read_rhos
    from maxlikelihood import bestparameters
    from plot_stats import pretty_rho
    import numpy as np

    a = b = n = 0; vara =  varb =  varn = 0; covab = covan = covbn = 0
    bestpar = bestparameters(samples)
    par_matcov = np.cov(samples) 
    if (par_matcov.size==1 ): variances = par_matcov
    else: variances = np.diagonal(par_matcov)
    covariances = sum( (par_matcov[i,i+1: ].tolist() for i in range(len(samples) - 1)) , [] )
    if(len(samples)==3):
        a, b, n = bestpar
        vara, varb, varn =  variances
        covab, covan, covbn =  covariances
    elif(len(samples)==2):
        a, b = bestpar
        vara, varb =  variances
        covab =  covariances[0]
    elif(len(samples)==1):
        a =  bestpar[0]
        vara =  variances
    else:
        print("Warning, test not defined")
    
    
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p,\
            sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4,\
            sig_rho5 = read_rhos(rhosfilename)

    #Ploting each term of the bias
    if(plots):
        xlim = [2., 300.]
        #supposing that a,b and n are idependent of rhos(scale independent)
        var0 = ((2*a*rho0p)**2)*vara +  (a**2)*(sig_rho0**2)
        var1 = ((2*b*rho1p)**2)*varb +  (b**2)*(sig_rho1**2)
        var2 = ((2*n*rho3p)**2)*varn +  (n**2)*(sig_rho3**2)
        varab =  vara*(b**2) + varb*(a**2) + 2*covab*(a*b)
        #varab = ((a*b)**2)*( (vara/((a)**2)) + (varb/((b)**2)) + 2*covab/(a*b) )
        var3 = 4*( (rho2p**2)*varab + (sig_rho2**2)*((a*b)**2)  )
        #var3 = 4*((a*b*rho2p)**2)*( varab/((a*b)**2) + (sig_rho2/rho2p)**2 )
        varbn =  varn*(b**2) + varb*(n**2) + 2*covbn*(b*n)
        #varbn = ((n*b)**2)*( (varn/((n)**2)) + (varb/((b)**2)) + 2*covbn/(b*n) ) 
        var4 = 4*( (rho4p**2)*varbn + (sig_rho4**2)*((n*b)**2)  )
        #var4 = 4*((n*b*rho4p)**2)*(varbn/((b*n)**2) + (sig_rho4/rho4p)**2)
        varan = varn*(a**2) + vara*(n**2) + 2*covan*(a*n)
        #varan = ((n*a)**2)*( (varn/((n)**2)) + (vara/((a)**2)) + 2*covan/(a*n) ) 
        var5 = 4*( (rho5p**2)*varan + (sig_rho5**2)*((n*a)**2)  )
        #var5 = 4*((n*a*rho5p)**2)*(varan/((a*n)**2) + (sig_rho5/rho5p)**2) 
        plt.clf()
        lfontsize = 7
        if (len(samples)==3):
            pretty_rho(meanr, (a**2)*rho0p, np.sqrt(var0), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b**2)*rho1p, np.sqrt(var1), legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (n**2)*rho3p, np.sqrt(var2), legend=r'$\eta^{2}\rho_{3}$', lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*a*b)*rho2p, np.sqrt(var3), legend=r'$2\alpha\beta \rho_{2}$',lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*b*n)*rho4p, np.sqrt(var4), legend=r'$2\beta\eta\rho_{4}$',lfontsize=lfontsize,  color='blue', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*n*a)*rho5p, np.sqrt(var5), legend=r'$2\eta\alpha\rho_{5}$', lfontsize=lfontsize, color='gray', ylabel='Correlations', xlim=xlim)
            print('Printing',  plotname)
            plt.savefig(plotname, dpi=200)
        if (len(samples)==2):
            pretty_rho(meanr, (a**2)*rho0p, np.sqrt(var0), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b**2)*rho1p, np.sqrt(var1), legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*a*b)*rho2p, np.sqrt(var3), legend=r'$2\alpha\beta \rho_{2}$',lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            print('Printing',  plotname)
            plt.savefig(plotname, dpi=200)
        if (len(samples)==1):
            pretty_rho(meanr, (a**2)*rho0p, np.sqrt(var0), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            print('Printing',  plotname)
            plt.savefig(plotname, dpi=200)
    
    #supposing that a,b and n are idependent of rhos(scale independent)
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
    return meanr, dxip, var_dxip
    
def RUNtest(args,  data, nwalkers, nsteps, i_guess, gflag, bflag,  eq='All',  moderr=False,  nsig=1,  namemc=None,  namecont=None, nameterms='terms_xip.png'):
    from chi2 import minimizeCHI2
    from maxlikelihood import MCMC, percentiles
    import numpy as np
    if (args.uwmprior):
        fitted_params =  np.array(i_guess)
        chisq = np.inf
    else:
        fitted_params, chisq = minimizeCHI2(data, i_guess, eq=eq,
                                            gflag=gflag,
                                            bflag=bflag,
                                            moderr=moderr)
    dof = len(data['rhos'][0])
    print("reduced Chi2:", chisq/dof )
    print("Found parameters" , fitted_params)
        
    samples = MCMC(fitted_params,data, nwalkers, nsteps, namemc,
                   namecont, eq=eq, gflag=gflag, bflag=bflag,
                   moderr=moderr, uwmprior=args.uwmprior,
                   plot=True )

    mcmcpars = percentiles(samples, nsig=nsig) 
    print('mcmc parameters',  mcmcpars)
    meanr, dxip, vardxip = getxipbias(samples, args.rhos, nameterms,  args.plots)
    return meanr, dxip, np.sqrt(vardxip) 

def main():
    import sys
    args = parse_args()
    sys.path.insert(0, args.srcpath)
    from plot_stats import plotallrhos,  plotalltaus,  pretty_rho
    from readjson import read_rhos, read_taus
    import numpy as np

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    plotspath = os.path.expanduser(args.outpath + 'plots/')
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(outpath): raise

        
    if (args.plots):
        xlim = [2., 300.]
        #Make directory where the ouput data will be 
        plotallrhos(args.rhos, outpath=plotspath, xlim=xlim)
        plotalltaus(args.taus, outpath=plotspath, xlim=xlim)
   
      
    # Save all the correlations in arrays
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p,\
    sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4,\
    sig_rho5 =read_rhos(args.rhos, args.maxscale)
    rhos = [rho0p, rho1p, rho2p, rho3p, rho4p, rho5p]
    sigrhos = [sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5]
    meanr2, tau0p, tau2p, tau5p, sig_tau0, sig_tau2, sig_tau5 = read_taus(args.taus, args.maxscale)
    taus = [tau0p, tau2p, tau5p]
    sigtaus = [sig_tau0, sig_tau2, sig_tau5]
    data = {}
    data['rhos'] = rhos
    data['sigrhos'] = sigrhos
    data['taus'] = taus
    data['sigtaus'] = sigtaus
    
    #Finding best alpha beta gamma
    nwalkers,  nsteps = 100,  10000
    moderr = False
    nsig = 1
    eq = 'All'
    i_guess0 = [ -0.01, 1, -1 ] #fiducial values

    if not (args.abe or args.ab or args.a): args.abe = True
    
    ## ALPHA-BETA-ETA
    if(args.abe):
        print("### Runing alpha, beta and eta test ### ")
        gflag, bflag = True, True
        i_guess = i_guess0
        namemc = plotspath + 'mcmc_alpha-beta-eta_eq_' + str(eq) + '_.png'
        namecont = plotspath +'contours_alpha-beta-eta_eq_' + str(eq) + '_.png'
        nameterms = plotspath +'termsdxip_alpha-beta-eta_eq_' + str(eq) + '_.png'
        vals = RUNtest(args, data, nwalkers, nsteps, i_guess, gflag, bflag, eq, moderr, nsig,  namemc, namecont, nameterms)
        writedxip( *vals,  filename=outpath + 'dxip-alpha-beta-eta.json' )
        writedxip_txt( *vals, outpath=outpath )
        if(args.plots):
            plt.clf()
            pretty_rho( *vals, legend=r"$\delta \xi_{+}$",  ylabel=r"$\delta \xi_{+}$",  xlim=xlim)
            fname = 'xibias_ABE.png' 
            print(plotspath +fname)
            plt.savefig(plotspath +fname, dpi=150)        
           
    ## ALPHA-BETA
    if(args.ab):
        print("### Runing alpha and beta test ### ")
        gflag, bflag = False, True
        i_guess = i_guess0[:2] #fiducial values
        namemc = plotspath + 'mcmc_alpha-beta_eq_' + str(eq) + '_.png'
        namecont = plotspath + 'contours_alpha-beta_eq_' + str(eq) + '_.png'
        nameterms = plotspath + 'termsdxip_alpha-beta_eq_' + str(eq) + '_.png'
        vals = RUNtest(args, data, nwalkers, nsteps, i_guess, gflag, bflag, eq, moderr, nsig,  namemc, namecont, nameterms)
        writedxip( *vals,  filename=outpath +'dxip-alpha-beta.json' )
        if(args.plots):
            plt.clf()
            pretty_rho( *vals, legend=r"$\delta \xi_{+}$",  ylabel=r"$\delta \xi_{+}$",  xlim=xlim)
            fname = 'xibias_AB.png' 
            print(plotspath+fname)
            plt.savefig(plotspath +fname, dpi=150)
               
    ## ALPHA
    if(args.a):
        print("### Runing alpha test ### ")
        gflag, bflag = False, False
        i_guess = i_guess0[:1] #fiducial values
        namemc = plotspath +'mcmc_alpha_eq_' + str(eq) + '_.png'
        namecont = plotspath +'contours_alpha_eq_' + str(eq) + '_.png'
        nameterms = plotspath +'termsdxip_alpha_eq_' + str(eq) + '_.png'
        vals = RUNtest(args, data, nwalkers, nsteps, i_guess, gflag, bflag, eq, moderr, nsig,  namemc, namecont, nameterms)
        writedxip( *vals,  filename=outpath+ 'dxip-alpha.json' )
        if(args.plots):
            plt.clf()
            pretty_rho( *vals, legend=r"$\delta \xi_{+}$",  ylabel=r"$\delta \xi_{+}$",  xlim=xlim)
            fname = 'xibias_A.png' 
            print(plotspath +fname)
            plt.savefig(plotspath +fname, dpi=150)
               
if __name__ == "__main__":
    main()
