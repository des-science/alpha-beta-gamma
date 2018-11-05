#Checking all the variances and different definitions.
#Ploting chisq by bin
def main():
    import os
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    from readjson import read_rhos, read_taus
    from chi2 import minimizeCHI2

    outpath = os.path.expanduser('/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/plots')
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

        
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 = read_rhos('/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_mod_epiff_magcut_irz.json')
    meanr2, tau0p, tau2p, tau5p, sig_tau0, sig_tau2, sig_tau5 =  read_taus('/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_irz.json')
    meanr3, tau0p_sn, tau2p_sn, tau5p_sn, sig_tau0_sn, sig_tau2_sn, sig_tau5_sn =  read_taus('/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_shapenoise_irz.json')

    plt.clf()
    plt.plot(meanr, sig_tau0**2, color='blue', label=r'$var(\tau_{0})$', marker='o')
    plt.plot(meanr, sig_tau2**2, color='red', label=r'$var(\tau_{2})$', marker='o')
    plt.plot(meanr, sig_tau5**2, color='green', label=r'$var(\tau_{5})$', marker='o')
    plt.plot(meanr, sig_tau0_sn**2, color='blue', label=r'$var(\tau_{0sn})$', marker='P')
    plt.plot(meanr, sig_tau2_sn**2, color='red', label=r'$var(\tau_{2sn})$', marker='P')
    plt.plot(meanr, sig_tau5_sn**2, color='green', label=r'$var(\tau_{5sn})$', marker='P')
    plt.legend(loc='upper right',  shadow=True,  fontsize=7)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    plt.ylim( [10**-25,10 **-9] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel('Variances', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()
    print("Printing :" + outpath +'/treecorr_vs_notes_variances.pdf')
    plt.savefig(outpath +'/treecorr_vs_notes_variances.pdf')

    plt.clf()
    plt.plot(meanr, (sig_tau0_sn**2 )/(sig_tau0**2), color='blue', label=r'$var(\tau_{0})/var(\tau_{0sn})$', marker='o')
    plt.plot(meanr, (sig_tau2_sn**2 )/(sig_tau2**2), color='green', label=r'$var(\tau_{2})/var(\tau_{2sn})$', marker='o')
    plt.plot(meanr, (sig_tau5_sn**2 )/(sig_tau5**2), color='red', label=r'$var(\tau_{5})/var(\tau_{5sn})$', marker='o')
    plt.legend(loc='best',  shadow=True,  fontsize=7)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    plt.ylim([0, 3])
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel('Variances', fontsize=24)
    plt.xscale('log')
    #plt.yscale('log', nonposy='clip')
    plt.tight_layout()
    print("Printing :" + outpath +'/ratio_treecorr_vs_notes_variances.pdf')
    plt.savefig(outpath +'/ratio_treecorr_vs_notes_variances.pdf')
    
    
    plt.clf()
    plt.plot(meanr, sig_tau0_sn**2, color='blue', label=r'$var(\tau_{0sn})$', marker='o')
    plt.plot(meanr, sig_tau2_sn**2, color='red', label=r'$var(\tau_{2sn})$', marker='o')
    plt.plot(meanr, sig_tau5_sn**2, color='green', label=r'$var(\tau_{5sn})$', marker='o')
    plt.plot(meanr, sig_rho0**2, color='black', label=r'$var(\rho_{0})$', marker='o')
    plt.plot(meanr, sig_rho1**2, color='yellow', label=r'$var(\rho_{1})$', marker='o')
    plt.plot(meanr, sig_rho2**2, color='gray', label=r'$var(\rho_{2})$', marker='o')
    plt.plot(meanr, sig_rho3**2, color='magenta', label=r'$var(\rho_{3})$', marker='o')
    plt.plot(meanr, sig_rho4**2, color='pink', label=r'$var(\rho_{4})$', marker='o')
    plt.plot(meanr, sig_rho5**2, color='orange', label=r'$var(\rho_{5})$', marker='o')
    plt.legend(loc='upper right',  shadow=True,  fontsize=7)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    plt.ylim( [10**-25,10 **-9] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel('Variances', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()
    print("Printing :" + outpath +'/all_variances_bybin.pdf' )
    plt.savefig(outpath +'/all_variances_bybin.pdf')
    













    
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
    
    
if __name__ == "__main__":
    main()
