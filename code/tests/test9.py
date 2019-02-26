#Plot the xipobs and the additional shift.
#plot each term of the ixpobs shift
#plot xi_teorico

def plotxipandbias(pars, pars_er, rhosfilename, tausfilename, xifilename, outpath):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    from readjson import read_rhos, read_taus,  read_xi
    from plot_stats import pretty_rho

    a, b, n = pars
    da, db, dn =  pars_er
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p,\
            sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4,\
            sig_rho5 = read_rhos(rhosfilename)
    meanr, tau0p, tau2p, tau5p, sig_tau0, sig_tau2,sig_tau5\
            = read_taus(tausfilename)
    
    meanr, xip_obs,  sig_xip = read_xi(xifilename)
    dxip = (a**2)*rho0p + (b**2)*rho1p + (n**2)*rho3p + (2*a*b)*rho2p + (2*b*n)*rho4p + (2*n*a)*rho5p
    sig_dxip = abs((a**2)*sig_rho0) + (b**2)*abs(sig_rho1) + (n**2)*abs(sig_rho3) + (2*a*b)*abs(sig_rho2) + (2*b*n)*abs(sig_rho4) + (2*n*a)*abs(sig_rho5)
    sig_dxip += 2*a*rho0p*da + 2*b*rho1p*db + 2*n*rho3p*dn + 2*rho2p*(a*db + da*b) + 2*rho4p*(b*dn + n*db) + 2*rho5p*(dn*a + n*da)                 
    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, xip_obs, sig_xip, sqrtn, legend=r'$\xi_{+}$',lfontsize=10,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, dxip, sig_dxip, sqrtn, legend=r'$\delta \xi_{+}$',lfontsize=10,  color='red', ylabel='Correlations', ylim=False)
    fname = 'xiobs_vs_xibias.png' 
    print(outpath +fname)
    plt.savefig(outpath +fname)

def plotbiasterms(pars, rhosfilename, tausfilename, outpath):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    from readjson import read_rhos, read_taus
    from plot_stats import pretty_rho

    a = pars[0];b = pars[1];n = pars[2]
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p,\
            sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4,\
            sig_rho5 = read_rhos(rhosfilename)
    meanr, tau0p, tau2p, tau5p, sig_tau0, sig_tau2,sig_tau5\
            = read_taus(tausfilename)

    sqrtn = 1
    plt.clf()
    pretty_rho(meanr, (a**2)*rho0p, sig_rho0, sqrtn, legend=r'$\alpha^{2} \rho_{0}$',lfontsize=10,  color='red', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, (b**2)*rho1p, sig_rho1, sqrtn, legend=r'$\beta^{2}\rho_{1}$',lfontsize=10,  color='green', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, (n**2)*rho3p, sig_rho3, sqrtn, legend=r'$\eta^{2}\rho_{3}$', lfontsize=10, color='black', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, (2*a*b)*rho2p, sig_rho2, sqrtn, legend=r'$2\alpha\beta \rho_{2}$',lfontsize=10,  color='yellow', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, (2*b*n)*rho4p, sig_rho4, sqrtn, legend=r'$2\beta\eta\rho_{4}$',lfontsize=10,  color='blue', ylabel='Correlations', ylim=False)
    pretty_rho(meanr, (2*n*a)*rho5p, sig_rho5, sqrtn, legend=r'$2\eta\alpha\rho_{5}$', lfontsize=10, color='gray', ylabel='Correlations', ylim=False)
    fname = 'xibias_parts.png' 
    print(outpath +fname)
    plt.savefig(outpath +fname)

def main():
    import os
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')

    outpath = "/home/dfa/sobreira/alsina/alpha-beta-gamma/code/tests/"
    
    rhosp =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_mod_epiff_magcut_sn_irz.json"
    tausp =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_irz.json"
    xip = "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/xi_xi_mod_riz.json"
    #theoretical xip
    theta = "/home/dfa/sobreira/alsina/alpha-beta-gamma/out/theta.txt"
    xip_teo = "/home/dfa/sobreira/alsina/alpha-beta-gamma/out/"

    alpha = 0.0183363720223933
    beta = 1.8480402601168726
    eta = -18.546803784998005
    da = 0.02868402123931059
    db = 1.5628513961084685
    de = 46.7407056515965
    pars = [alpha, beta, eta]
    pars_er = [da, db, de]
    
    
    plotbiasterms(pars, rhosp, tausp, outpath)
    plotxipandbias(pars, pars_er, rhosp, tausp, xip, outpath)
if __name__ == "__main__":
    main()
