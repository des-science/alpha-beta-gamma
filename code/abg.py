import os

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Alpha beta gamma test solving the system of equatiosn and plotting correlations')
    
    parser.add_argument('--rsgcorr',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_irz.json',
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
    sigtaus = [tau5p, sig_tau0, sig_tau2, sig_tau5]
    data = {}
    data['rhos'] = rhos
    data['sigrhos'] = sigrhos
    data['taus'] = taus
    data['sigtaus'] = sigtaus

    eq = 0
    
    dof = len(rhos[0])
    ## ALPHA
    gflag, bflag = False, False
    i_guess = [0] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    print("alpha" , fitted_params)
    print("Chi2 reduced:", chisq/dof )
    x_arr= -0.05, 0.05, 100
    filename1 = outpath+'/chisq_only_alpha.pdf'
    filename2 = outpath+'/chisqshifted_only_alpha.pdf'
    plotCHI2(None, data,x_arr,filename1, eq=eq,gflag=gflag,bflag=bflag)
    plotCHI2shifted(None, data, x_arr, chisq, filename2, eq=eq, gflag=gflag, bflag=bflag)
    #MCMC(fitted_params,data, svalue=chisq, eq=eq, gflag=gflag, bflag=bflag )
    #OneParMaxLike(fitted_param,data, svalue=svalue,eq=eq, gflag=gflag,bflag = bflag)
    
    ## ALPHA-BETA
    gflag, bflag = False, True
    i_guess = [0,-1] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    print("alpha, beta" , fitted_params)
    print("Chi2 reduced:", chisq/dof )
    namemc = outpath+'/mcmc_alpha-beta.pdf'
    namecont = outpath+'/contours_alpha-beta.pdf'
    MCMC(fitted_params,data ,namemc, namecont, eq=eq, gflag=gflag, bflag=bflag )

    ## ALPHA-BETA-GAMMA
    gflag, bflag = True, True
    i_guess = [0,-1,- 1] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    print("alpha, beta, gamma:" , fitted_params)
    print("Chi2 reduced:", chisq/dof )
    namemc = outpath+'/mcmc_alpha-beta-eta.pdf'
    namecont = outpath+'/contours_alpha-beta-eta.pdf'
    MCMC(fitted_params,data, namemc, namecont,  eq=eq, gflag=gflag, bflag=bflag )

    
if __name__ == "__main__":
    main()
