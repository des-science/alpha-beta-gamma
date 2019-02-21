#run abn test and save files arrays of parameters going out to a maximum bin. 
#plotting each term in the equation of correlations, to see if there are reason to cancel some part of the model.
def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/tests/',
                        help='location of the output of the files')

    args = parser.parse_args()

    return args
def write_pars(name, meanr, a_c, a_l, a_r, b_c, b_l, b_r, d_c, d_l, d_r):
    import json
    stats = [
        meanr.tolist(), 
        a_c,
        a_l,
        a_r,
        b_c,
        b_l,
        b_r,
        d_c,
        d_l,
        d_r,
    ]
    print('stat_file = ', name)
    with open(name,'w') as fp:
        json.dump([stats], fp)
    print('Done writing ',name)

def read_pars(name):
    import json
    import numpy as np
    with open(name,'r') as f:
        stats = json.load(f)
    if len(stats) == 1:
        stats = stats[0]
        ( meanr, 
          a_c,
          a_l,
          a_r,
          b_c,
          b_l,
          b_r,
          d_c,
          d_l,
          d_r,
        )=stats[:10]
    a_c = np.array(a_c); a_l = np.array(a_l); a_r = np.array(a_r)
    a_c = np.array(a_c); a_l = np.array(a_l); a_r = np.array(a_r)
    a_c = np.array(a_c); a_l = np.array(a_l); a_r = np.array(a_r)
    return meanr, a_c, a_l, a_r, b_c, b_l, b_r, d_c, d_l,d_r
    
def run_parspatch(outpath):
    from readjson import read_rhos, read_taus
    from chi2 import minimizeCHI2
    from maxlikelihood import MCMC, percentiles
    import numpy as np

    nwalkers,  nsteps = 100,  1000
    eq = None; moderr = False
    gflag, bflag = True, True
    i_guess = [0,-1,- 1] #fiducial values
    data = {}
    for patch in range(1, 5):
        rhosp =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_mod_patch" + str(patch) + "_irz.json"
        tausp =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tau_all_galaxy-reserved_mod_patch" + str(patch) + "_irz.json"

        a_c = []; a_l = []; a_r = []
        b_c = []; b_l = []; b_r = []
        d_c = []; d_l = []; d_r = []
        for ibin in range(1, 33):
            meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p,\
            sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4,\
            sig_rho5 = read_rhos(rhosp, maxbin=ibin)
       
            meanr, tau0p, tau2p, tau5p, sig_tau0, sig_tau2,sig_tau5 =\
            read_taus(tausp, maxbin=ibin)
  
            rhos = [rho0p, rho1p, rho2p, rho3p, rho4p, rho5p]
            sigrhos = [sig_rho0, sig_rho1, sig_rho2, sig_rho3,
                       sig_rho4, sig_rho5]
            taus = [tau0p, tau2p, tau5p]
            sigtaus = [sig_tau0, sig_tau2, sig_tau5]
            data['rhos'] = rhos
            data['sigrhos'] = sigrhos
            data['taus'] = taus
            data['sigtaus'] = sigtaus

            fit_pars, chisq = minimizeCHI2(data, i_guess, eq=eq,
                                           gflag=gflag, bflag=bflag,
                                           moderr=moderr)
            samples = MCMC(fit_pars,data, nwalkers, nsteps, eq=eq,
                            gflag=gflag, bflag=bflag, moderr=moderr,
                            plot=False )
            mcmcpars = percentiles(samples, nsig=2) 
            a_c.append(mcmcpars[0][0]); a_l.append(mcmcpars[0][1]);a_r.append(mcmcpars[0][2])
            b_c.append(mcmcpars[1][0]); b_l.append(mcmcpars[1][1]);b_r.append(mcmcpars[1][2])
            d_c.append(mcmcpars[2][0]); d_l.append(mcmcpars[2][1]);d_r.append(mcmcpars[2][2])
        write_pars(outpath + 'parspatch' + str(patch) + '.json', meanr, a_c, a_l, a_r, b_c, b_l, b_r, d_c, d_l, d_r)

def plotlog(outpath):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    colors = ['black', 'green', 'blue', 'red', 'gray', 'pink']
    for patch in range(1, 5):
        meanr, a_c, a_l, a_r, b_c, b_l, b_r, d_c, d_l, d_r = read_pars(outpath + 'parspatch' + str(patch) + '.json')
        
        label = "P" + str(patch)
        width = np.diff(meanr)*0.25
        width = np.append(width, width[-1:])
        dx = meanr + (patch - 1)*width
        plt.figure(0)
        colormat=np.where(np.array(a_c)>0, 'r','y')
        plt.bar(dx, np.absolute(a_c), yerr=[ -np.array(a_l),
                                             np.array(a_r)],
                width=width,
                color=colors[patch],
                label=label, log=True,
                ec="k", align="edge",
                capsize=0)
        #plt.scatter(dx+ width*0.5, np.absolute(a_c), marker="D",
        #            alpha=0.8, color=colormat, zorder=2)
        
        #plt.plot(meanr, a_c, color=colors[patch],  label=label, marker='o', linestyle='None')
        #plt.plot(meanr, -np.array(a_c), color=colors[patch],  marker='v', linestyle='None')
        #plt.errorbar(meanr, a_c , yerr=[ -np.array(a_l),  np.array(a_r)] , color=colors[patch], capsize=5)
        #plt.errorbar(meanr, -np.array(a_c) ,yerr=[ -np.array(a_l),  np.array(a_r)] , color=colors[patch], capsize=5) 
        plt.figure(1)
        colormat=np.where(np.array(b_c)>0, 'r','y')
        plt.bar(dx, np.absolute(b_c), yerr=[ -np.array(b_l),
                                             np.array(b_r)],
                width=width,
                color=colors[patch],
                label=label, log=True,
                ec="k", align="edge",
                capsize=5)
        plt.scatter(dx+ width*0.5, np.absolute(b_c), marker="D",
                    alpha=0.8, color=colormat, zorder=2)
        #plt.plot(meanr, b_c, color=colors[patch], label=label, marker='o')
        #plt.plot(meanr, -np.array(b_c), color=colors[patch], marker='v')
        #plt.errorbar(meanr, b_c , yerr=[ -np.array(b_l),  np.array(b_r)] , color=colors[patch], capsize=10)
        #plt.errorbar(meanr, -np.array(b_c) , yerr=[-np.array(b_l),np.array(b_r)] , color=colors[patch], capsize=10)
        plt.figure(2)
        colormat=np.where(np.array(d_c)>0, 'r','y')
        plt.bar(dx, np.absolute(d_c), yerr=[ -np.array(d_l),
                                             np.array(d_r)],
                width=width,
                color=colors[patch],
                label=label, log=True,
                ec="k", align="edge",
                capsize=5)
        plt.scatter(dx+ width*0.5, np.absolute(d_c), marker="D",
                    alpha=0.8, color=colormat, zorder=2)
        #plt.plot(meanr, d_c, color=colors[patch], label=label, marker='o')
        #plt.plot(meanr, -np.array(d_c), color=colors[patch], marker='v')
        #plt.errorbar(meanr, d_c , yerr=[ -np.array(d_l),  np.array(d_r)] , color=colors[patch],  capsize=10)
        #plt.errorbar(meanr, -np.array(d_c) , yerr=[-np.array(d_l), np.array(d_r)] , color=colors[patch],  capsize=10)
    plt.figure(0)
    plt.legend(loc='best', fontsize=10)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\alpha$', fontsize=24)
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim( [1.e-5, 1.e+1] )
    plt.tight_layout()
    print("Printing :", outpath +'alpha_quadrants.pdf')
    plt.savefig(outpath +'alpha_quadrants.pdf')
    print("Printed :", outpath +'alpha_quadrants.pdf')
    plt.figure(1)
    plt.legend(loc='best', fontsize=10)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\beta$', fontsize=24)
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim( [1.e-3, 1.e+3] )
    plt.tight_layout()
    print("Printing :", outpath +'/beta_quadrants.pdf')
    plt.savefig(outpath +'/beta_quadrants.pdf')
    print("Printed :", outpath +'/beta_quadrants.pdf')
    plt.figure(2)
    plt.legend(loc='best', fontsize=10)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\eta$', fontsize=24)
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim( [1.e-3, 1.e+3] )
    plt.tight_layout()
    print("Printing :", outpath +'/eta_quadrants.pdf')
    plt.savefig(outpath +'/eta_quadrants.pdf')
    print("Printed :", outpath +'/eta_quadrants.pdf')
        
def plotlineal(outpath,  bar=False):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    colors = ['black', 'green', 'blue', 'red', 'gray', 'pink']
    for patch in range(1, 5):
        meanr, a_c, a_l, a_r, b_c, b_l, b_r, d_c, d_l, d_r = read_pars(outpath + 'parspatch' + str(patch) + '.json')
        
        label = "P" + str(patch)
        width = np.diff(meanr)*0.25
        width = np.append(width, width[-1:])
        dx = meanr + (patch - 1)*width
        plt.figure(0)
        if(bar):
            plt.bar(dx, a_c, yerr=[ -np.array(a_l), np.array(a_r)],
                    width=width, color=colors[patch], label=label,
                    ec="k", align="edge", capsize=0)
        else:
            plt.scatter(meanr, a_c, color=colors[patch],  label=label, marker='o', s=10)
            plt.errorbar(meanr, a_c , yerr=[ -np.array(a_l),
                                             np.array(a_r)] ,
                         color=colors[patch],
                         capsize=0, linestyle='')
        
        plt.figure(1)
        if(bar):
            plt.bar(dx, b_c, yerr=[ -np.array(b_l), np.array(b_r)],
                    width=width, color=colors[patch], label=label,
                    ec="k", align="edge", capsize=0)
        else:
            plt.scatter(meanr, b_c, color=colors[patch],  label=label, marker='o', s=10)
            plt.errorbar(meanr, b_c , yerr=[ -np.array(b_l),
                                             np.array(b_r)] ,
                         color=colors[patch],
                         capsize=0, linestyle='')
        
        plt.figure(2)
        if(bar):
            plt.bar(dx, d_c, yerr=[ -np.array(d_l), np.array(d_r)],
                    width=width, color=colors[patch], label=label,
                    ec="k", align="edge", capsize=0)
        else:
            plt.scatter(meanr, d_c, color=colors[patch],  label=label, marker='o', s=10)
            plt.errorbar(meanr, d_c , yerr=[ -np.array(d_l),
                                             np.array(d_r)] ,
                         color=colors[patch],
                         capsize=0, linestyle='')
        
    plt.figure(0)
    plt.legend(loc='best', fontsize=10)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\alpha$', fontsize=24)
    plt.xscale('log')
    plt.xlim( [ 4e-1, 300] )
    plt.ylim( [ - 0.4, 0.4] )
    plt.tight_layout()
    print("Printing :", outpath +'alpha_quadrants.png')
    plt.savefig(outpath +'alpha_quadrants.png')
    print("Printed :", outpath +'alpha_quadrants.png')
    plt.figure(1)
    plt.legend(loc='lower right', fontsize=10)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\beta$', fontsize=24)
    plt.xscale('log')
    plt.xlim( [ 4e-1, 300] )
    plt.ylim( [ - 75, 40] )
    plt.tight_layout()
    print("Printing :", outpath +'beta_quadrants.png')
    plt.savefig(outpath +'/beta_quadrants.png')
    print("Printed :", outpath +'beta_quadrants.png')
    plt.figure(2)
    plt.legend(loc='best', fontsize=10)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\eta$', fontsize=24)
    plt.xscale('log')
    plt.xlim( [ 4e-1, 300] )
    plt.ylim( [ -600, 600] )
    plt.tight_layout()
    print("Printing :", outpath +'eta_quadrants.png')
    plt.savefig(outpath +'eta_quadrants.png')
    print("Printed :", outpath +'eta_quadrants.png')
        
def main():
    import os
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    
    
    args = parse_args()
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise
    
    run_parspatch(outpath)
    #plotlog(outpath)
    plotlineal(outpath, bar=False)
    


    
if __name__ == "__main__":
    main()
