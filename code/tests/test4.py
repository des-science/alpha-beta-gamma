#Comparing rhostats modified and unmodified
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



    #Comparing mod and unmod tau correlations
    modtaus =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_irz.json"
    unmodtaus =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_unmod_irz.json"
    meanr, tau0p , tau2p , tau5p , sig_tau0 , sig_tau2 , sig_tau5  = read_taus(modtaus)
    meanr_m, tau0p_m, tau2p_m, tau5p_m, sig_tau0_m, sig_tau2_m, sig_tau5_m = read_taus(unmodtaus)
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr,tau0p ,  sig_tau0 , sqrtn, legend=r'$\tau_{0}$', color='blue', ylim=False, lfontsize=10)
    pretty_rho(meanr_m, tau0p_m,  sig_tau0_m, sqrtn, legend=r'$\tau_0*$', color='blue', marker='P', ylim=False, lfontsize=10)
    pretty_rho(meanr,tau2p ,  sig_tau2 , sqrtn, legend=r'$\tau_{2}$', color='green', ylim=False, lfontsize=10)
    pretty_rho(meanr_m, tau2p_m,  sig_tau2_m, sqrtn,  legend=r'$\tau_2*$', color='green', marker='P', ylim=False, lfontsize=10)
    pretty_rho(meanr,tau5p ,  sig_tau5 , sqrtn, legend=r'$\tau_{5}$', color='red', ylim=False, lfontsize=10)
    pretty_rho(meanr_m, tau5p_m,  sig_tau5_m, sqrtn, legend=r'$\tau_5*$', color='red', marker='P', ylim=False, lfontsize=10, ylabel=r'$\tau(\theta)$')
    plt.xlim([0, 1000])
    print("Printing :", outpath +'/taustatsmeanvsnomean1.pdf')
    plt.savefig(outpath +'/taustatsmeanvsnomean1.pdf')
    

    #Comparing mod and unmod rhos correlations
    modrhoseobs =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_mod_eobs_magcut_irz.json"
    modrhosepiff =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_unmod_eobs_magcut_irz.json"
    meanr_obs, rho0p_obs, rho1p_obs, rho2p_obs, rho3p_obs, rho4p_obs, rho5p_obs, sig_rho0_obs, sig_rho1_obs, sig_rho2_obs, sig_rho3_obs, sig_rho4_obs, sig_rho5_obs = read_rhos(modrhoseobs)
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 = read_rhos(modrhosepiff)
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr,rho1p_obs,  sig_rho1_obs, sqrtn, legend=r'$\rho_{1}$', color='blue', ylim=False, lfontsize=10)
    pretty_rho(meanr, rho1p,  sig_rho1, sqrtn, legend=r'$\rho_1*$', color='blue', marker='P', ylim=False, lfontsize=10)
    pretty_rho(meanr,rho3p_obs,  sig_rho3_obs, sqrtn, legend=r'$\rho_{3}$', color='green', ylim=False, lfontsize=10)
    pretty_rho(meanr, rho3p,  sig_rho3, sqrtn,  legend=r'$\rho_3*$', color='green', marker='P', ylim=False, lfontsize=10)
    pretty_rho(meanr,rho4p_obs,  sig_rho4_obs, sqrtn, legend=r'$\rho_{4}$', color='red', ylim=False, lfontsize=10)
    pretty_rho(meanr, rho4p,  sig_rho4, sqrtn, legend=r'$\rho_4*$', color='red', marker='P', ylim=False, lfontsize=10)
    plt.xlim([0, 400])
    print("Printing :", outpath +'/statsmeanvsnomean1.pdf')
    plt.savefig(outpath +'/statsmeanvsnomean1.pdf')
    plt.clf()
    pretty_rho(meanr,rho2p_obs,  sig_rho2_obs, sqrtn, legend=r'$\rho_{2}$', color='blue', ylim=False, lfontsize=10)
    pretty_rho(meanr, rho2p,  sig_rho2, sqrtn, legend=r'$\rho_2*$', color='blue', marker='P', ylim=False, lfontsize=10)
    pretty_rho(meanr,rho5p_obs,  sig_rho5_obs, sqrtn, legend=r'$\rho_{5}$', color='green', ylim=False, lfontsize=10)
    pretty_rho(meanr, rho5p,  sig_rho5, sqrtn,  legend=r'$\rho_5*$', color='green', marker='P', ylim=False, lfontsize=10)
    plt.xlim([0, 400])
    print("Printing :", outpath +'/statsmeanvsnomean2.pdf')
    plt.savefig(outpath +'/statsmeanvsnomean2.pdf')
    plt.clf()
    pretty_rho(meanr,rho0p_obs,  sig_rho0_obs, sqrtn, legend=r'$\rho_{0}$', color='blue', ylim=False, lfontsize=10)
    pretty_rho(meanr, rho0p,  sig_rho0, sqrtn, legend=r'$\rho_0*$', color='blue', marker='P', ylim=False, lfontsize=10)
    plt.xlim([0, 400])
    print("Printing :", outpath +'/statsmeanvsnomean0.pdf')
    plt.savefig(outpath +'/statsmeanvsnomean0.pdf')
        
    #Reading a ploting reserved stars correlations
    modrhoseobs =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_mod_eobs_magcut_irz.json"
    modrhosepiff =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_mod_epiff_magcut_irz.json"
    meanr_obs, rho0p_obs, rho1p_obs, rho2p_obs, rho3p_obs, rho4p_obs, rho5p_obs, sig_rho0_obs, sig_rho1_obs, sig_rho2_obs, sig_rho3_obs, sig_rho4_obs, sig_rho5_obs = read_rhos(modrhoseobs)
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 = read_rhos(modrhosepiff)
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr,rho1p_obs,  sig_rho1_obs, sqrtn, legend=r'$\rho_{1}(e_{obs})$', color='blue')
    pretty_rho(meanr, rho1p,  sig_rho1, sqrtn, legend=r'$\rho_1(e_{piff})$', color='green', marker='P')
    print("Printing :", outpath +'/deltarho1_modeobs.pdf')
    plt.savefig(outpath +'/deltarho1_modeobs.pdf')


    plt.clf()
    pretty_rho(meanr, rho2p + rho1p, np.sqrt(sig_rho1**2 + sig_rho2**2), sqrtn, legend=r'$\rho_2(e_{piff})+\rho_1(e_{piff})$', color='blue')
    pretty_rho(meanr,rho2p_obs, sig_rho2_obs, sqrtn, legend=r'$\rho_2(e_{obs})$', color='green', marker='P')
    print("Printing :", outpath +'/r2+r1modepiff.pdf')
    plt.savefig(outpath +'/r2+r1modepiff.pdf')

if __name__ == "__main__":
    main()
   
