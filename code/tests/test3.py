#plotting each term in the equation of correlations, to see if there are reason to cancel some part of the model.
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='')
    
    parser.add_argument('--rsgcorr',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_irz.json',
                        help='Json file with the reserved stars - galaxies correlations')
    parser.add_argument('--rsrscorr',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_mod_epiff_magcut_irz.json',
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
    from plot_stats import pretty_rho
    from chi2 import minimizeCHI2
    import numpy as np

    args = parse_args()
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

        
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 = read_rhos(args.rsrscorr)
    meanr2, tau0p, tau2p, tau5p, sig_tau0, sig_tau2, sig_tau5 =  read_taus(args.rsgcorr)
    
    dof = len(meanr)
    rhos = [rho0p, rho1p, rho2p, rho3p, rho4p, rho5p]
    sigrhos = [sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5]
    taus = [tau0p, tau2p, tau5p]
    sigtaus = [sig_tau0, sig_tau2, sig_tau5]
    data = {}
    data['rhos'] = rhos
    data['sigrhos'] = sigrhos
    data['taus'] = taus
    data['sigtaus'] = sigtaus

    ### EQUATION 0 ###
    #ALPHA-BETA-GAMMA
    eq = 0
    gflag, bflag = True, True
    i_guess = [0,-1,- 1] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    alpha, beta, eta = fitted_params
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, tau0p, sig_tau0, sqrtn, legend=r'$\tau_{0}$', lfontsize=15,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, alpha*rho0p, sig_rho0, sqrtn, legend=r'$\alpha \rho_{0}$',lfontsize=15,  color='red', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, beta*rho2p, sig_rho2, sqrtn, legend=r'$\beta\rho_{2}$',lfontsize=15,  color='green', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, eta*rho5p, sig_rho5, sqrtn, legend=r'$\eta\rho_{5}$', lfontsize=15, color='black', ylabel='Correlations', ylim=False)
    plt.title(r'$\chi^{2}_{\nu}$ =' + str(chisq / dof)[:6], fontsize=5)
    print(outpath +'/eq0_abn_corr_pars0.pdf')
    plt.savefig(outpath +'/eq0_abn_corr_pars0.pdf')
    #ALPHA-BETA
    gflag, bflag = False, True
    i_guess = [0,-1] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    alpha, beta = fitted_params
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, tau0p, sig_tau0, sqrtn, legend=r'$\tau_{0}$', lfontsize=15,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, alpha*rho0p, sig_rho0, sqrtn, legend=r'$\alpha \rho_{0}$',lfontsize=15,  color='red', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, beta*rho2p, sig_rho2, sqrtn, legend=r'$\beta\rho_{2}$',lfontsize=15,  color='green', ylabel='Correlations', ylim=False)
    plt.title(r'$\chi^{2}_{\nu}$ =' + str(chisq / dof)[:6], fontsize=5)
    print(outpath +'/eq0_ab_corr_pars0.pdf')
    plt.savefig(outpath +'/eq0_ab_corr_pars0.pdf')
    #ALPHA-ONLY
    gflag, bflag = False, False
    i_guess = [0] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    alpha = fitted_params
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, tau0p, sig_tau0, sqrtn, legend=r'$\tau_{0}$', lfontsize=15,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, alpha*rho0p, sig_rho0, sqrtn, legend=r'$\alpha \rho_{0}$',lfontsize=15,  color='red', ylabel='Correlations', ylim=False)
    plt.title(r'$\chi^{2}_{\nu}$ =' + str(chisq / dof)[:6], fontsize=5)
    print(outpath +'/eq0_a_corr_pars0.pdf')
    plt.savefig(outpath +'/eq0_a_corr_pars0.pdf')


    ### EQUATION 1 ###
    #ALPHA-BETA-GAMMA
    eq = 1
    gflag, bflag = True, True
    i_guess = [0,-1,- 1] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    alpha, beta, eta = fitted_params
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, tau2p, sig_tau2, sqrtn, legend=r'$\tau_{2}$', lfontsize=15,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, alpha*rho2p, sig_rho2, sqrtn, legend=r'$\alpha \rho_{2}$',lfontsize=15,  color='red', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, beta*rho1p, sig_rho1, sqrtn, legend=r'$\beta\rho_{1}$',lfontsize=15,  color='green', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, eta*rho4p, sig_rho4, sqrtn, legend=r'$\eta\rho_{4}$', lfontsize=15, color='black', ylabel='Correlations', ylim=False)
    plt.title(r'$\chi^{2}_{\nu}$ =' + str(chisq / dof)[:6], fontsize=5)
    print(outpath +'/eq1_abn_corr_pars0.pdf')
    plt.savefig(outpath +'/eq1_abn_corr_pars0.pdf')
    #ALPHA-BETA
    gflag, bflag = False, True
    i_guess = [0,-1] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    alpha, beta = fitted_params
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, tau2p, sig_tau2, sqrtn, legend=r'$\tau_{2}$', lfontsize=15,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, alpha*rho2p, sig_rho2, sqrtn, legend=r'$\alpha \rho_{2}$',lfontsize=15,  color='red', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, beta*rho1p, sig_rho1, sqrtn, legend=r'$\beta\rho_{1}$',lfontsize=15,  color='green', ylabel='Correlations', ylim=False)
    plt.title(r'$\chi^{2}_{\nu}$ =' + str(chisq / dof)[:6], fontsize=5)
    print(outpath +'/eq1_ab_corr_pars0.pdf')
    plt.savefig(outpath +'/eq1_ab_corr_pars0.pdf')
    #ALPHA-ONLY
    gflag, bflag = False, False
    i_guess = [0] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    alpha = fitted_params
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, tau2p, sig_tau2, sqrtn, legend=r'$\tau_{2}$', lfontsize=15,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, alpha*rho2p, sig_rho2, sqrtn, legend=r'$\alpha \rho_{2}$',lfontsize=15,  color='red', ylabel='Correlations', ylim=False)
    plt.title(r'$\chi^{2}_{\nu}$ =' + str(chisq / dof)[:6], fontsize=5)
    print(outpath +'/eq1_a_corr_pars0.pdf')
    plt.savefig(outpath +'/eq1_a_corr_pars0.pdf')


    ### EQUATION 2 ###
    #ALPHA-BETA-GAMMA
    eq = 2
    gflag, bflag = True, True
    i_guess = [0,-1,- 1] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    alpha, beta, eta = fitted_params
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, tau5p, sig_tau5, sqrtn, legend=r'$\tau_{5}$', lfontsize=15,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, alpha*rho5p, sig_rho5, sqrtn, legend=r'$\alpha \rho_{5}$',lfontsize=15,  color='red', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, beta*rho4p, sig_rho4, sqrtn, legend=r'$\beta\rho_{4}$',lfontsize=15,  color='green', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, eta*rho3p, sig_rho3, sqrtn, legend=r'$\eta\rho_{3}$', lfontsize=15, color='black', ylabel='Correlations', ylim=False)
    plt.title(r'$\chi^{2}_{\nu}$ =' + str(chisq / dof)[:6], fontsize=5)
    print(outpath +'/eq2_abn_corr_pars0.pdf')
    plt.savefig(outpath +'/eq2_abn_corr_pars0.pdf')
    #ALPHA-BETA
    gflag, bflag = False, True
    i_guess = [0,-1] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    alpha, beta = fitted_params
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, tau5p, sig_tau5, sqrtn, legend=r'$\tau_{5}$', lfontsize=15,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, alpha*rho5p, sig_rho5, sqrtn, legend=r'$\alpha \rho_{5}$',lfontsize=15,  color='red', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, beta*rho4p, sig_rho4, sqrtn, legend=r'$\beta\rho_{4}$',lfontsize=15,  color='green', ylabel='Correlations', ylim=False)
    plt.title(r'$\chi^{2}_{\nu}$ =' + str(chisq / dof)[:6], fontsize=5)
    print(outpath +'/eq2_ab_corr_pars0.pdf')
    plt.savefig(outpath +'/eq2_ab_corr_pars0.pdf')
    #ALPHA-ONLY
    gflag, bflag = False, False
    i_guess = [0] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    alpha = fitted_params
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, tau5p, sig_tau5, sqrtn, legend=r'$\tau_{5}$', lfontsize=15,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, alpha*rho5p, sig_rho5, sqrtn, legend=r'$\alpha \rho_{5}$',lfontsize=15,  color='red', ylabel='Correlations', ylim=False)
    plt.title(r'$\chi^{2}_{\nu}$ =' + str(chisq / dof)[:6], fontsize=5)
    print(outpath +'/eq2_a_corr_pars0.pdf')
    plt.savefig(outpath +'/eq2_a_corr_pars0.pdf')
    
    

    ### All system of EQUATIONS ###
    eq = None
    #ALPHA-BETA-GAMMA
    gflag, bflag = True, True
    i_guess = [0,-1,- 1] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
  
    alpha, beta, eta = fitted_params
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, tau0p + tau2p + tau5p, np.sqrt(sig_tau0**2 + sig_tau2**2 + sig_tau5**2) , sqrtn, legend=r'$\tau_{0}+\tau_{2}+\tau_{5}$', lfontsize=15,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, alpha*(rho0p+rho2p+rho5p) , sig_rho0, sqrtn, legend=r'$\alpha\times(\rho_{0}+\rho_{2}+\rho_{5})$',lfontsize=15,  color='red', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, beta*(rho2p +rho1p + rho4p), sig_rho2, sqrtn, legend=r'$\beta\times(\rho_{2}+\rho_{1}+\rho_{4})$',lfontsize=15,  color='green', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, eta*(rho5p+rho4p + rho3p) , sig_rho5, sqrtn, legend=r'$\eta\times(\rho_{5}+\rho_{4}+\rho_{3})$', lfontsize=15, color='black', ylabel='Correlations', ylim=False)
    plt.title(r'$\chi^{2}_{\nu}$ =' + str(chisq / dof)[:6], fontsize=5)
    print(outpath +'/soe_abg_corr_pars0.pdf')
    plt.savefig(outpath +'/soe_abg_corr_pars0.pdf')

    #ALPHA-BETA
    gflag, bflag = False, True
    i_guess = [0,-1] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    alpha, beta = fitted_params
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, tau0p + tau2p + tau5p, np.sqrt(sig_tau0**2 + sig_tau2**2 + sig_tau5**2) , sqrtn, legend=r'$\tau_{0}+\tau_{2}+\tau_{5}$', lfontsize=15,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, alpha*(rho0p+rho2p+rho5p) , sig_rho0, sqrtn, legend=r'$\alpha\times(\rho_{0}+\rho_{2}+\rho_{5})$',lfontsize=15,  color='red', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, beta*(rho2p +rho1p + rho4p), sig_rho2, sqrtn, legend=r'$\beta\times(\rho_{2}+\rho_{1}+\rho_{4})$',lfontsize=15,  color='green', ylabel='Correlations', ylim=False)
    plt.title(r'$\chi^{2}_{\nu}$ =' + str(chisq / dof)[:6], fontsize=5)
    print(outpath +'/soe_ab_corr_pars0.pdf')
    plt.savefig(outpath +'/soe_ab_corr_pars0.pdf')

    #ALPHA
    gflag, bflag = False, False
    i_guess = [0] #fiducial values
    fitted_params, chisq =  minimizeCHI2(data, i_guess,  eq=eq, gflag=gflag, bflag=bflag)
    alpha = fitted_params
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, tau0p + tau2p + tau5p, np.sqrt(sig_tau0**2 + sig_tau2**2 + sig_tau5**2) , sqrtn, legend=r'$\tau_{0}+\tau_{2}+\tau_{5}$', lfontsize=15,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, alpha*(rho0p+rho2p+rho5p) , sig_rho0, sqrtn, legend=r'$\alpha\times(\rho_{0}+\rho_{2}+\rho_{5})$',lfontsize=15,  color='red', ylabel='Correlations', ylim=False)
    plt.title(r'$\chi^{2}_{\nu}$ =' + str(chisq / dof)[:6], fontsize=5)
    print(outpath +'/soe_a_corr_pars0.pdf')
    plt.savefig(outpath +'/soe_a_corr_pars0.pdf')
    

    
    
    
if __name__ == "__main__":
    main()
