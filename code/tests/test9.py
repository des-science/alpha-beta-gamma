#Create fit file with all the biases and the covariance matrix of dxip
# not implelement cross redshift bins, nor scale bins.
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--rhos',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_irz_26-03-19_all_reserved_mod_epiff_magcut_sn.json',
                        help='Json file with the reserved stars -reserved stars correlations')
    parser.add_argument('--tausfolder',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tomo_taus/',help='location of the folder containing all the taus files. They must endin the number of the patch')
    parser.add_argument('--outpath', default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/', help='location of the output of the files')
    parser.add_argument('--filename', default='abg_dxip_tomo.fits', help='Name of the fit file where info of dxip will be saved ')
    parser.add_argument('--maxscale', default=15,  type=float, 
                        help='Limit the analysis to certain maximum scale, units are determined by .json file with the correlations')
    parser.add_argument('--uwmprior', default=False,
                        action='store_const', const=True, help='Run all tomographic correlations')
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
def findbinfile(ipath, i, j, symm=False):
    import numpy as np
    files = np.array(os.listdir(ipath))
    filesnoext= [os.path.splitext(file)[0] for file in files  ] 
    b1 = np.array([ f.endswith(str(i)+'_'+str(j)) for f in filesnoext ])
    b2 = np.array([ f.endswith(str(j)+'_'+str(i)) for f in filesnoext ])
    if(symm):
        out_file = files[(b1 | b2)]
        if (len(out_file)>1):
            print('WARNING: bin file ', str(i), str(j) , 'is repeated, probably you 2pt is not simmetric')
        return (out_file[0])
    else:
        out_file = files[b1]
        if (len(out_file)!=1):
            print('WARNING: bin file ', str(i), str(j) , 'is repeated or does not exist')
            return None
        else:
            return (out_file[0])
  
def getxipbias(samples, rhosfilename, plotname=None,  plots=False):
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
  
def RUNtest(args,  data, nwalkers, nsteps, i_guess, gflag, bflag,  eq='All',  moderr=False,  nsig=1,  namemc=None,  namecont=None, nameterms=None):
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
                   plot=False )

    mcmcpars = percentiles(samples, nsig=nsig) 
    print('mcmc parameters',  mcmcpars)
    meanr, dxip, vardxip = getxipbias(samples, args.rhos, nameterms)
    return meanr, dxip, np.sqrt(vardxip) 

def write_fit(data, names, filename):
    import fitsio
    from fitsio import FITS,FITSHDR
    import os.path
    if not os.path.isfile(filename):
        fits = FITS(filename,'rw')
        fits.write(data, names=names, clobber=False)
    else:
        fits = FITS(filename,'rw')
        fits[-1].append(data)
def main():
    import sys
    args = parse_args()
    sys.path.insert(0, args.srcpath)
    import numpy as np
    import itertools
    from readjson import read_rhos, read_taus
    import fitsio
    from fitsio import FITS,FITSHDR
    from astropy.io import fits
    
    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
        

         
    nwalkers,  nsteps = 10,  1000
    moderr = False
    nsig = 1
    eq = 'All'
    i_guess0 = [ -0.01, 1, -1 ] #fiducial values

    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p,\
    sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4,\
    sig_rho5 =read_rhos(args.rhos, args.maxscale)
    rhos = [rho0p, rho1p, rho2p, rho3p, rho4p, rho5p]
    sigrhos = [sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5]
    data = {}
    data['rhos'] = rhos
    data['sigrhos'] = sigrhos

    names=['BIN1','BIN2','ANGBIN', 'ANG', 'VALUE']
    forms = ['i4', 'i4', 'i4',  'f8',  'f8']
    dtype = dict(names = names, formats=forms)
    nrows = 20
    outdata = np.recarray((nrows, ), dtype=dtype)

    nbins = 4
    a=[i for i in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.combinations_with_replacement(a, 2):
        bin_pairs.append(p)


    covmat = np.zeros(shape=(nbins*nrows , nbins*nrows ))
 
    covmatbin = None; listofmat = []
    for i,j in bin_pairs:            
        taufilename = findbinfile(args.tausfolder, i, j)
        if taufilename is None:
            continue
        print(i, j, taufilename)

        meanr2, tau0p, tau2p, tau5p, sig_tau0, sig_tau2, sig_tau5 \
            = read_taus(args.tausfolder + taufilename, args.maxscale)
        taus = [tau0p, tau2p, tau5p]
        sigtaus = [sig_tau0, sig_tau2, sig_tau5]
        
        data['taus'] = taus
        data['sigtaus'] = sigtaus
        
        if not (args.abe or args.ab or args.a): args.abe = True   
    
        ## ALPHA-BETA-ETA
        if(args.abe):
            print("### Runing alpha, beta and eta test ### ")
            gflag, bflag = True, True
            i_guess = i_guess0
            meanr, dxip, vardxip = RUNtest(args, data, nwalkers,
                                           nsteps, i_guess, gflag, bflag, eq, moderr,
                                           nsig)

            covmat[nrows*(i-1):nrows*i,nrows*(j-1) :nrows*j] = np.diag(vardxip)
            
            angarray = meanr
            valuearray =  dxip
            bin1array = np.array([i]*len(meanr))
            bin2array = np.array([j]*len(meanr))
            angbinarray = np.array([i for i in range(len(meanr))])
            array_list = [bin1array, bin2array, angbinarray, angarray,  valuearray ]
            for array, name in zip(array_list, names): outdata[name] = array 
            write_fit(outdata, names, args.filename)
        
            
        ## ALPHA-BETA
        if(args.ab):
            print("### Runing alpha and beta test ### ")
            gflag, bflag = False, True
            i_guess = i_guess0[:2] #fiducial values
            meanr, dxip, vardxip = RUNtest(args, data, nwalkers,
                                           nsteps, i_guess, gflag, bflag, eq, moderr,
                                           nsig)

            
        ## ALPHA
        if(args.a):
            print("### Runing alpha test ### ")
            gflag, bflag = False, False
            i_guess = i_guess0[:1] #fiducial values
            meanr, dxip, vardxip = RUNtest(args, data, nwalkers,
                                           nsteps, i_guess, gflag, bflag, eq, moderr,
                                           nsig) 

    hdulist = fits.open(args.filename)
    hdulist[1].name = 'xip'
    covmathdu = fits.ImageHDU(covmat, name='COVMAT')
    hdulist.insert(1, covmathdu) 
    hdulist.writeto(args.filename, clobber=True)
    print(covmat)
    
if __name__ == "__main__":
    main()
    
#hstacks = [ np.block(i) for i in [listofmat[j*nblocks:(j+1)*nblocks] for j in range(nblocks)] ]
#covmat =  np.vstack(hstacks)
#for i in hstacks: print(i.shape)

