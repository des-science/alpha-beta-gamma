import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Alpha beta gamma test solving the fitting problem of system ofequatiosn, plotting correlations and final correlation function withbias')
    
    parser.add_argument('--taus', nargs='+',
                        default=['/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/TAU0P_bin_0.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/TAU2P_bin_0.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/TAU5P_bin_0.fits'],
                        help='Fits file containing all rho stats fitsfiles. --taus tau0.fits tau2.fits tau5.fits')
    parser.add_argument('--rhos', nargs='+',
                        default=['/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO0P.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO1P.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO2P.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO3P.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO4P.fits',
                                 '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO5P.fits'],
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
    parser.add_argument('--eq', default=4, type=int, 
                        help='Select equations to be used for istance --eq=0, 4 represent the whole system of equations')
    parser.add_argument('--abe', default=False,
                        action='store_const', const=True, help='Run alpha, beta and eta test.')
    parser.add_argument('--ab', default=False,
                        action='store_const', const=True, help='Run alpha and beta test.')
    parser.add_argument('--ae', default=False,
                        action='store_const', const=True, help='Run alpha and eta test.')
    parser.add_argument('--be', default=False,
                        action='store_const', const=True, help='Run beta and eta test.')
    parser.add_argument('--a', default=False,
                        action='store_const', const=True, help='Run only alpha test.')
    parser.add_argument('--b', default=False,
                        action='store_const', const=True, help='Run only beta test.')
    parser.add_argument('--e', default=False,
                        action='store_const', const=True, help='Run only eta test.')
    
    
    
    args = parser.parse_args()

    return args
def corrmatrix(cov):
    import numpy as np
    cov = np.mat(cov)
    D = np.diag(np.sqrt(np.diag(cov)))
    d = np.linalg.inv(D)
    corr = d*cov*d
    return corr
def writexipbias(samples,rhonames, nsig=1, plots=False,xim=False,nameterms='terms_dxi.png',dxiname='dxi.png',namecovmat='covm_pars.png',filename='dxi.fits'):
    from readjson import read_rhos
    from totalmaxlikelihood import bestparameters, percentiles
    from plot_stats import pretty_rho
    from readfits import read_corr
    from astropy.io import fits
    import numpy as np

    mcmcpars = percentiles(samples, nsig=nsig) 
    print('nsig=', nsig, ' mcmc parameters: ',  mcmcpars)
    
    ##Format of the fit file output
    names=['BIN1', 'BIN2','ANGBIN', 'VALUE', 'ANG']
    forms = ['i4', 'i4', 'i4',  'f8',  'f8']
    dtype = dict(names = names, formats=forms)
    nrows = 20
    outdata = np.recarray((nrows, ), dtype=dtype)
    namesout=['TAU0P', 'TAU2P', 'TAU5P', 'TAU0M','TAU2M', 'TAU5M']
    
    #plot covariance matrix of parameters alpha, beta and eta.
    if plots:
        par_matcov = np.cov(samples)
        corr=corrmatrix(par_matcov)
        #print(par_matcov)
        #print(corr)
        cov_vmin=np.min(corr)
        plt.imshow(corr,cmap='viridis'+'_r', interpolation='nearest',
                   aspect='auto', origin='lower', vmin=cov_vmin, vmax=1.)
        plt.colorbar()
        plt.title(r'$\alpha \mid \beta \mid \eta $')
        plt.savefig(namecovmat, dpi=500)
        print(namecovmat, 'Printed!')

    a = b = n = 0; vara =  varb =  varn = 0; covab = covan = covbn = 0
    bestpar = bestparameters(samples)
    print("Best pars", bestpar)
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
        var0 = ((2*a*rho0)**2)*vara +  (a**2)*(sig_rho0**2)
        var1 = ((2*b*rho1)**2)*varb +  (b**2)*(sig_rho1**2)
        var2 = ((2*n*rho3)**2)*varn +  (n**2)*(sig_rho3**2)
        varab =  vara*(b**2) + varb*(a**2) + 2*covab*(a*b)
        #varab = ((a*b)**2)*( (vara/((a)**2)) + (varb/((b)**2)) + 2*covab/(a*b) )
        var3 = 4*( (rho2**2)*varab + (sig_rho2**2)*((a*b)**2)  )
        #var3 = 4*((a*b*rho2p)**2)*( varab/((a*b)**2) + (sig_rho2/rho2p)**2 )
        varbn =  varn*(b**2) + varb*(n**2) + 2*covbn*(b*n)
        #varbn = ((n*b)**2)*( (varn/((n)**2)) + (varb/((b)**2)) + 2*covbn/(b*n) ) 
        var4 = 4*( (rho4**2)*varbn + (sig_rho4**2)*((n*b)**2)  )
        #var4 = 4*((n*b*rho4p)**2)*(varbn/((b*n)**2) + (sig_rho4/rho4p)**2)
        varan = varn*(a**2) + vara*(n**2) + 2*covan*(a*n)
        #varan = ((n*a)**2)*( (varn/((n)**2)) + (vara/((a)**2)) + 2*covan/(a*n) ) 
        var5 = 4*( (rho5**2)*varan + (sig_rho5**2)*((n*a)**2)  )
        #var5 = 4*((n*a*rho5p)**2)*(varan/((a*n)**2) + (sig_rho5/rho5p)**2) 
        plt.clf()
        lfontsize = 7
        if (len(samples)==3):
            pretty_rho(meanr, (a**2)*rho0, np.sqrt(np.diag(cov_rho0)), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b**2)*rho1, np.sqrt(var1), legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (n**2)*rho3, np.sqrt(var2), legend=r'$\eta^{2}\rho_{3}$', lfontsize=lfontsize, color='black', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*a*b)*rho2, np.sqrt(var3), legend=r'$2\alpha\beta \rho_{2}$',lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*b*n)*rho4, np.sqrt(var4), legend=r'$2\beta\eta\rho_{4}$',lfontsize=lfontsize,  color='blue', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*n*a)*rho5, np.sqrt(var5), legend=r'$2\eta\alpha\rho_{5}$', lfontsize=lfontsize, color='gray', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (len(samples)==2):
            pretty_rho(meanr, (a**2)*rho0, np.sqrt(var0), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (b**2)*rho1, np.sqrt(var1), legend=r'$\beta^{2}\rho_{1}$',lfontsize=lfontsize,  color='green', ylabel='Correlations', xlim=xlim)
            pretty_rho(meanr, (2*a*b)*rho2, np.sqrt(var3), legend=r'$2\alpha\beta \rho_{2}$',lfontsize=lfontsize,  color='yellow', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
        if (len(samples)==1):
            pretty_rho(meanr, (a**2)*rho0, np.sqrt(var0), legend=r'$\alpha^{2} \rho_{0}$',lfontsize=lfontsize,  color='red', ylabel='Correlations', xlim=xlim)
            print('Printing',  nameterms)
            plt.savefig(nameterms, dpi=200)
    
    #supposing that a,b and n are idependent of rhos(scale independent)
    dxi = (a**2)*rho0 + (b**2)*rho1 + (n**2)*rho3 + (2*a*b)*rho2 + (2*b*n)*rho4 + (2*n*a)*rho5
    f1 = 2*(a*rho0 + b*rho2 + n*rho5)     
    f2 = 2*(b*rho1 + a*rho2 + n*rho4)
    f3 = 2*(n*rho3 + b*rho4 + a*rho5)
    f4 = a**2 ; f5 = b**2; f6 = 2*a*b
    f7 = n**2 ; f8 = 2*b*n; f9 = 2*n*a 
    covmat_dxi = np.diag( (f1**2)*vara + (f2**2)*varb + (f3**2)*varn + + 2*(f1*f2*covab + f1*f3*covan + f2*f3*covbn) ) \
    + (f4**2)*(cov_rho0) + (f5**2)*(cov_rho1) + (f6**2)*(cov_rho2) + (f7**2)*(cov_rho3) +(f8**2)*(cov_rho4) + (f9**2)*(cov_rho5) 

    if(plots):
        plt.clf()
        pretty_rho(meanr, dxi, np.sqrt(np.diag(covmat_dxi)) , legend=r"$\delta \xi_{+}$",  ylabel=r"$\delta \xi_{+}$",  xlim=xlim)
        print('Printing',  dxiname)
        plt.savefig(dxiname, dpi=150)

    nrows = len(dxi)
    hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu])
    covmathdu = fits.ImageHDU(covmat_dxi, name='COVMAT')
    hdul.insert(1, covmathdu)
    angarray = meanr
    valuearray =  np.array(dxi)
    bin1array = np.array([ -999]*nrows)
    bin2array = np.array([ -999]*nrows)
    angbinarray = np.arange(nrows)
    array_list = [bin1array, bin2array, angbinarray, valuearray,  angarray ]
    for array, name in zip(array_list, names): outdata[name] = array 
    corrhdu = fits.BinTableHDU(outdata, name='xi')
    hdul.insert(2, corrhdu)
    hdul.writeto(filename, clobber=True)
    print(filename,'Written!')

def RUNtest(i_guess, data, nwalkers, nsteps,  eq='All', mflags=[True, True, True] ,   moderr=False, uwmprior=False, minimize=True ):
    from totalchi2 import minimizeCHI2
    from totalmaxlikelihood import MCMC
    import numpy as np
    if (uwmprior or (not minimize)):
        iguess =  np.array(i_guess)
        chisq = np.inf
    if(minimize):
        fitted_params, chisq = minimizeCHI2(data, i_guess, eq=eq,
                                            mflags=mflags,
                                            moderr=moderr)
        dof = len(data['rhos'][0])
        print("reduced Chi2:", chisq/dof )
        print("Found parameters" , fitted_params)
        iguess = fitted_params
    
    samples,  chains = MCMC(i_guess,data, nwalkers, nsteps, eq=eq,
                   mflags=mflags, moderr=moderr, uwmprior=uwmprior)
   
    return samples,  chains
    
def main():
    import sys
    args = parse_args()
    sys.path.insert(0, args.srcpath)
    from readfits import read_corr
    from plot_stats import plotallrhosfits,  plotallrhoscorrmatfits, plotalltausfits,  plotalltauscorrmatfits,  plot_samplesdist
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
        ['/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO0M.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO1M.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO2M.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO3M.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO4M.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/RHO5M.fits']
        args.taus = \
        ['/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/TAU0M.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/TAU2M.fits',
         '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/TAU5M.fits']

    if (args.plots):
        xlim = [2., 300.]
        ylims = [[1.e-11,1.e-7 ],[3.e-8 ,3.e-4 ]]
        rhostitle = ''
        plotallrhosfits(args.rhos, outpath=plotspath, title=rhostitle, xlim=xlim, ylims=ylims)
        plotallrhoscorrmatfits(args.rhos, outpath=plotspath)
        plotalltausfits(args.taus, outpath=plotspath, xlim=xlim)
        plotalltauscorrmatfits(args.taus, outpath=plotspath)

    rhonames = args.rhos
    meanr, rho0, cov_rho0 = read_corr(rhonames[0], maxscale=args.maxscale)
    meanr, rho1, cov_rho1 = read_corr(rhonames[1], maxscale=args.maxscale)
    meanr, rho2, cov_rho2 = read_corr(rhonames[2], maxscale=args.maxscale)
    meanr, rho3, cov_rho3 = read_corr(rhonames[3], maxscale=args.maxscale)
    meanr, rho4, cov_rho4 = read_corr(rhonames[4], maxscale=args.maxscale)
    meanr, rho5, cov_rho5 = read_corr(rhonames[5], maxscale=args.maxscale)

    rhos = [rho0, rho1, rho2, rho3, rho4, rho5]
    covrhos = [cov_rho0, cov_rho1, cov_rho2, cov_rho3, cov_rho4, cov_rho5]
    taunames = args.taus
    meanr, tau0, cov_tau0 = read_corr(taunames[0], maxscale=args.maxscale)
    meanr, tau2, cov_tau2 = read_corr(taunames[1], maxscale=args.maxscale)
    meanr, tau5, cov_tau5 = read_corr(taunames[2], maxscale=args.maxscale)
    taus = [tau0, tau2, tau5]
    covtaus = [cov_tau0, cov_tau2, cov_tau5]
    data = {}
    data['rhos'] = rhos
    data['covrhos'] = covrhos
    data['taus'] = taus
    data['covtaus'] = covtaus
    
    #for i in taus: print i.size
    #for i in covrhos: print i.shape

    #Finding best alpha beta gamma
    nwalkers,  nsteps = 100,  1000
    moderr = False
    minimize = True
    nsig = 1
    eq = args.eq
    print("Using equations: ", eq)
    i_guess0 = [ -0.01, 1, -1 ] #fiducial values

    if not (args.abe or args.ab or args.a): args.abe = True
    
    ## ALPHA-BETA-ETA
    if(args.abe):
        print("### Runing alpha, beta and eta test ### ")
        mflags = [True, True, True] ##alpha,beta,eta
        namemc = plotspath + 'mcmc_alpha-beta-eta_eq_' + str(eq) + '_.png'
        namecont = plotspath +'contours_alpha-beta-eta_eq_' + str(eq) + '_.png'
        nameterms = plotspath +'termsdxip_alpha-beta-eta_eq_' + str(eq) + '_.png'
        namecovmat = plotspath +'covmatrix_alpha-beta-eta_eq_' + str(eq) + '_.png'
        namedxip = plotspath +'xibias_abe_' + str(eq) + '_.png'
        filename =  outpath +'abe_dxi.fits'            
    ## ALPHA-BETA
    if(args.ab):
        print("### Runing alpha and beta test ### ")
        mflags = [True, True, False] ##alpha,beta,eta
        namemc = plotspath + 'mcmc_alpha-beta_eq_' + str(eq) + '_.png'
        namecont = plotspath + 'contours_alpha-beta_eq_' + str(eq) + '_.png'
        nameterms = plotspath + 'termsdxip_alpha-beta_eq_' + str(eq) + '_.png'
        namecovmat = plotspath +'covmatrix_alpha-beta_eq_' + str(eq) + '_.png'
        namedxip = plotspath +'xibias_ab_' + str(eq) + '_.png'
        filename =  outpath +'ab_dxi.fits'
    ## ALPHA-ETA
    if(args.ae):
        print("### Runing alpha and eta test ### ")
        mflags = [True, False, True] ##alpha,eta,eta
        namemc = plotspath + 'mcmc_alpha-eta_eq_' + str(eq) + '_.png'
        namecont = plotspath + 'contours_alpha-eta_eq_' + str(eq) + '_.png'
        nameterms = plotspath + 'termsdxip_alpha-eta_eq_' + str(eq) + '_.png'
        namecovmat = plotspath +'covmatrix_alpha-eta_eq_' + str(eq) + '_.png'
        namedxip = plotspath +'xibias_ae_' + str(eq) + '_.png'
        filename =  outpath +'ae_dxi.fits'
    ## BETA-ETA
    if(args.be):
        print("### Runing beta and eta test ### ")
        mflags = [True, False, True] ##beta,eta,eta
        namemc = plotspath + 'mcmc_beta-eta_eq_' + str(eq) + '_.png'
        namecont = plotspath + 'contours_beta-eta_eq_' + str(eq) + '_.png'
        nameterms = plotspath + 'termsdxip_beta-eta_eq_' + str(eq) + '_.png'
        namecovmat = plotspath +'covmatrix_beta-eta_eq_' + str(eq) + '_.png'
        namedxip = plotspath +'xibias_be_' + str(eq) + '_.png'
        filename =  outpath +'be_dxi.fits' 
    ## ALPHA
    if(args.a):
        print("### Runing alpha test ### ")
        mflags = [True, False, False] ##alpha,beta,eta
        namemc = plotspath +'mcmc_alpha_eq_' + str(eq) + '_.png'
        namecont = plotspath +'contours_alpha_eq_' + str(eq) + '_.png'
        nameterms = plotspath +'termsdxip_alpha_eq_' + str(eq) + '_.png'
        namecovmat = plotspath +'covmatrix_alpha_eq_' + str(eq) + '_.png'
        namedxip = plotspath +'xibias_a_' + str(eq) + '_.png'
        filename =  outpath +'a_dxi.fits'
    ## Beta
    if(args.b):
        print("### Runing beta test ### ")
        mflags = [False, True, False] ##alpha,beta,eta
        namemc = plotspath +'mcmc_beta_eq_' + str(eq) + '_.png'
        namecont = plotspath +'contours_beta_eq_' + str(eq) + '_.png'
        nameterms = plotspath +'termsdxip_beta_eq_' + str(eq) + '_.png'
        namecovmat = plotspath +'covmatrix_beta_eq_' + str(eq) + '_.png'
        namedxip = plotspath +'xibias_b_' + str(eq) + '_.png'
        filename =  outpath +'b_dxi.fits'
    ## Eta
    if(args.e):
        print("### Runing eta test ### ")
        mflags = [False, False, True] ##alpha,eta,eta
        namemc = plotspath +'mcmc_eta_eq_' + str(eq) + '_.png'
        namecont = plotspath +'contours_eta_eq_' + str(eq) + '_.png'
        nameterms = plotspath +'termsdxip_eta_eq_' + str(eq) + '_.png'
        namecovmat = plotspath +'covmatrix_eta_eq_' + str(eq) + '_.png'
        namedxip = plotspath +'xibias_e_' + str(eq) + '_.png'
        filename =  outpath +'e_dxi.fits'

    i_guess = np.array(i_guess0)[np.array(mflags)].tolist()
    samples, chains = RUNtest(i_guess, data, nwalkers, nsteps, eq=eq,
                      mflags=mflags, moderr=moderr,
                      uwmprior=args.uwmprior, minimize= minimize)
    #samples= np.c_[[par[int(0.2 * len(par)):] for par in samples]].T
    #print("Total samples", [len(i) for i in samples] )
    if(args.plots): plot_samplesdist(samples, chains, mflags, nwalkers, nsteps,  namemc, namecont )
    writexipbias(samples, args.rhos, plots=args.plots, xim=args.xim,
                 nameterms=nameterms, dxiname=namedxip,
                 namecovmat=namecovmat, filename=filename )
    
  
    
if __name__ == "__main__":
    main()

