import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Alpha beta gamma test solving the fitting problem of system ofequatiosn, plotting correlations and final correlation function withbias')
    
    parser.add_argument('--taus', nargs='+',
                        default=['/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/TAU0p.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/TAU2p.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/TAU5p.fits'],
                        help='Fits file containing all rho stats fitsfiles. --taus tau0.fits tau2.fits tau5.fits')
    parser.add_argument('--rhos', nargs='+',
                        default=['/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO0p.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO1p.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO2p.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO3p.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO4p.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO5p.fits'],
                        help='Fits file containing all rho stats fitsfiles. --rhos rho0.fits rho1.fits rho2.fits.. rho5.fits')
    parser.add_argument('--xim', default=False,
                        action='store_const', const=True,
                        help='trecorr return xim instead of xip')
    parser.add_argument('--maxscale', default=None, type=float, 
                        help='Limit the analysis to certain maximum scale, units are determined by .json file with the correlations')
    parser.add_argument('--uwmprior', default=False,
                        action='store_const', const=True, help='Run all tomographic correlations')
    parser.add_argument('--plots', default=False,
                        action='store_const', const=True, help='Plot correlations functions')
    parser.add_argument('--outpath',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/',
                        help='location of the output of the final contaminant')
    parser.add_argument('--plotspath',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/plots/',
                        help='location of the plots.')
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
def corrmatrix(cov):
    import numpy as np
    cov = np.mat(cov)
    D = np.diag(np.sqrt(np.diag(cov)))
    d = np.linalg.inv(D)
    corr = d*cov*d
    return corr
def writexipbias(samples,rhosfilenames,xim=False,plots=False,nameterms='terms_dxi.png',namedxip='dxi.png',namecovmat='covm_pars.png',filename='dxip.fits'):
    from readjson import read_rhos
    from maxlikelihood import bestparameters
    from plot_stats import pretty_rho
    from readfits import read_corr
    from astropy.io import fits
    import numpy as np

    #plot covariance matrix of parameters alpha, beta and eta.
    if plots:
        par_matcov = np.cov(samples)
        corr=corrmatrix(par_matcov)
        print(par_matcov)
        print(corr)
        cov_vmin=np.min(corr)
        plt.imshow(corr,cmap='viridis'+'_r', interpolation='nearest',
                   aspect='auto', origin='lower', vmin=cov_vmin, vmax=1.)
        plt.colorbar()
        plt.title(r'$\alpha \mid \beta \mid \eta $')
        plt.savefig(namecovmat, dpi=500)
        print(namecovmat, 'Printed!')

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
        print("Warning, test type not defined")
    
    
    rhonames = args.rhos
    meanr, rho0, cov_rho0 = read_corr(rhonames[0])
    meanr, rho1, cov_rho1 = read_corr(rhonames[1])
    meanr, rho2, cov_rho2 = read_corr(rhonames[2])
    meanr, rho3, cov_rho3 = read_corr(rhonames[3])
    meanr, rho4, cov_rho4 = read_corr(rhonames[4])
    meanr, rho5, cov_rho5 = read_corr(rhonames[5])
    sig_rho0 =  np.sqrt(np.diag(cov_rho0))
    sig_rho1 =  np.sqrt(np.diag(cov_rho1))
    sig_rho2 =  np.sqrt(np.diag(cov_rho2))
    sig_rho3 =  np.sqrt(np.diag(cov_rho3))
    sig_rho4 =  np.sqrt(np.diag(cov_rho4))
    sig_rho5 =  np.sqrt(np.diag(cov_rho5))

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
            pretty_rho(meanr, (a**2)*rho0p, np.sqrt(np.diag(cov_rho0)), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b**2)*rho1p, np.sqrt(var1), legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (n**2)*rho3p, np.sqrt(var2), legend=r'$\eta^{2}\rho_{3}$', lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*a*b)*rho2p, np.sqrt(var3), legend=r'$2\alpha\beta \rho_{2}$',lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*b*n)*rho4p, np.sqrt(var4), legend=r'$2\beta\eta\rho_{4}$',lfontsize=lfontsize,  color='blue', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*n*a)*rho5p, np.sqrt(var5), legend=r'$2\eta\alpha\rho_{5}$', lfontsize=lfontsize, color='gray', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (len(samples)==2):
            pretty_rho(meanr, (a**2)*rho0p, np.sqrt(var0), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b**2)*rho1p, np.sqrt(var1), legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*a*b)*rho2p, np.sqrt(var3), legend=r'$2\alpha\beta \rho_{2}$',lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (len(samples)==1):
            pretty_rho(meanr, (a**2)*rho0p, np.sqrt(var0), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
    
    #supposing that a,b and n are idependent of rhos(scale independent)
    dxip = (a**2)*rho0p + (b**2)*rho1p + (n**2)*rho3p + (2*a*b)*rho2p + (2*b*n)*rho4p + (2*n*a)*rho5p
    f1 = 2*(a*rho0p + b*rho2p + n*rho5p)     
    f2 = 2*(b*rho1p + a*rho2p + n*rho4p)
    f3 = 2*(n*rho3p + b*rho4p + a*rho5p)
    f4 = a**2 ; f5 = b**2; f6 = 2*a*b
    f7 = n**2 ; f8 = 2*b*n; f9 = 2*n*a 
    covmat_dxip = np.diag( (f1**2)*vara + (f2**2)*varb + (f3**2)*varn + + 2*(f1*f2*covab + f1*f3*covan + f2*f3*covbn) ) \
    + (f4**2)*(cov_rho0) + (f5**2)*(cov_rho1) + (f6**2)*(cov_rho2) + (f7**2)*(cov_rho3) +(f8**2)*(cov_rho4) + (f9**2)*(cov_rho5) 

    if(plots):
        plt.clf()
        pretty_rho(meanr, dxip, np.sqrt(np.diag(covmat_dxip)) , legend=r"$\delta \xi_{+}$",  ylabel=r"$\delta \xi_{+}$",  xlim=xlim)
        print('Printing',  dxipname)
        plt.savefig(dxipname, dpi=150)

    nrows = len(dxip)
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])
    covmathdu = fits.ImageHDU(covmat_dxip, name='COVMAT')
    hdul.insert(1, covmathdu)
    angarray = meanr
    valuearray =  np.array(dxip)
    bin1array = np.array([ -999]*nrows)
    bin2array = np.array([ -999]*nrows)
    angbinarray = np.arange(nrows)
    array_list = [bin1array, bin2array, angbinarray, valuearray,  angarray ]
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name=nam)
    hdul.insert(2, corrhdu)
    if xim:
        hdul.writeto(filename + 'm.fits', clobber=True)
    else:
        hdul.writeto(filename + 'p.fits', clobber=True)
     
def RUNtest(args,  data, nwalkers, nsteps, i_guess, gflag, bflag,  eq='All',  moderr=False,  nsig=1,  namemc=None,  namecont=None):
    from fullchi2 import minimizeCHI2
    from fullmaxlikelihood import MCMC, percentiles
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
    return samples
    
def main():
    import sys
    args = parse_args()
    sys.path.insert(0, args.srcpath)
    from readfits import read_corr
    from plot_stats import plotallrhosfits,  plotallrhoscorrmatfits, plotalltausfits,  plotalltauscorrmatfits
    import numpy as np

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    plotspath = os.path.expanduser(args.plotspath)
    try:
        if not os.path.exists(plotspath):
            os.makedirs(plotspath)
    except OSError:
        if not os.path.exists(outpath): raise

    if(args.xim):
        args.rhos = \
        ['/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO0m.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO1m.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO2m.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO3m.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO4m.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO5m.fits']
        args.taus = \
        ['/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/TAU0m.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/TAU2m.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/TAU5m.fits']

    if (args.plots):
        xlim = [2., 300.]
        ylims = [[4.e-5, 3.e-4],[1.e-11,1.e-7 ],[3.e-8 ,1.e-6 ]]
        rhostitle = 'Jackknife errors'
        plotallrhosfits(args.rhos, outpath=plotspath, title=rhostitle, xlim=xlim, ylims=ylims)
        plotallrhoscorrmatfits(args.rhos, outpath=plotspath)
        #plotalltausfits(args.taus, outpath=plotspath, xlim=xlim)
        #plotalltauscorrmatfits(args.taus, outpath=plotspath)

    rhonames = args.rhos
    meanr, rho0, cov_rho0 = read_corr(rhonames[0])
    meanr, rho1, cov_rho1 = read_corr(rhonames[1])
    meanr, rho2, cov_rho2 = read_corr(rhonames[2])
    meanr, rho3, cov_rho3 = read_corr(rhonames[3])
    meanr, rho4, cov_rho4 = read_corr(rhonames[4])
    meanr, rho5, cov_rho5 = read_corr(rhonames[5])
    rhos = [rho0, rho1, rho2, rho3, rho4, rho5]
    covrhos = [cov_rho0, cov_rho1, cov_rho2, cov_rho3, cov_rho4, cov_rho5]

    taunames = args.tau
    meanr, tau0, cov_tau0 = read_corr(taunames[0])
    meanr, tau2, cov_tau2 = read_corr(taunames[1])
    meanr, tau5, cov_tau5 = read_corr(taunames[2])
    taus = [tau0, tau2, tau5]
    covtaus = [cov_tau0, cov_tau2, cov_tau5]
    data = {}
    data['rhos'] = rhos
    data['covrhos'] = covrhos
    data['taus'] = taus
    data['covtaus'] = covtaus

    #Finding best alpha beta gamma
    nwalkers,  nsteps = 100,  1000
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
        namecovmat = plotspath +'covmatrix_alpha-beta-eta_eq_' + str(eq) + '_.png'
        namedxip = plotspath +'xobias_abe_' + str(eq) + '_.png'
        filename =  outspath +'abe_dxip.fits'
        samples = RUNtest(args, data, nwalkers, nsteps, i_guess, gflag, bflag, eq, moderr, nsig,  namemc, namecont)
        writexipbias(samples, args.rhos, args.plots, args.xim,  nameterms, namedxip, namecovmat, filename )
            
    ## ALPHA-BETA
    if(args.ab):
        print("### Runing alpha and beta test ### ")
        gflag, bflag = False, True
        i_guess = i_guess0[:2] #fiducial values
        namemc = plotspath + 'mcmc_alpha-beta_eq_' + str(eq) + '_.png'
        namecont = plotspath + 'contours_alpha-beta_eq_' + str(eq) + '_.png'
        nameterms = plotspath + 'termsdxip_alpha-beta_eq_' + str(eq) + '_.png'
        namecovmat = plotspath +'covmatrix_alpha-beta_eq_' + str(eq) + '_.png'
        namedxip = plotspath +'xibias_ab_' + str(eq) + '_.png'
        filename =  outspath +'ab_dxip.fits'
        samples = RUNtest(args, data, nwalkers, nsteps, i_guess, gflag, bflag, eq, moderr, nsig,  namemc, namecont)
        writexipbias(samples, args.rhos, args.plots, args.xim,  nameterms, namedxip, namecovmat, filename )
       
               
    ## ALPHA
    if(args.a):
        print("### Runing alpha test ### ")
        gflag, bflag = False, False
        i_guess = i_guess0[:1] #fiducial values
        namemc = plotspath +'mcmc_alpha_eq_' + str(eq) + '_.png'
        namecont = plotspath +'contours_alpha_eq_' + str(eq) + '_.png'
        nameterms = plotspath +'termsdxip_alpha_eq_' + str(eq) + '_.png'
        namecovmat = plotspath +'covmatrix_alpha-beta-eta_eq_' + str(eq) + '_.png'
        namedxip = plotspath +'xibias_a_' + str(eq) + '_.png'
        filename =  outspath +'a_dxip.fits'
        samples = RUNtest(args, data, nwalkers, nsteps, i_guess, gflag, bflag, eq, moderr, nsig,  namemc, namecont)
        writexipbias(samples, args.rhos, args.plots, args.xim,  nameterms, namedxip, namecovmat, filename )
    
if __name__ == "__main__":
    main()

