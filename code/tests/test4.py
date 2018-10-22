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
   
