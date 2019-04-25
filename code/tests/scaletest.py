#run abn test and save files arrays of parameters going out to a maximum bin. 
#plotting each term in the equation of correlations, to see if there are reason to cancel some part of the model.
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--rhos',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/allrhos_4jk.fits',
                        help='location of the file containing all the rhos by patch.')
    parser.add_argument('--taus',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/correlations/alltaus_4jk.fits',
                        help='location of the file containing all the taus by patch')
    parser.add_argument('--outpath', default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/plots/',
                        help='location of the output of the files, i.e plots.')

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
    
def findpatchfile(ipath, i):
    import numpy as np
    files = np.array(os.listdir(ipath))
    filesnoext= [os.path.splitext(file)[0] for file in files  ] 
    b1 = np.array([ f.endswith(str(i)) for f in filesnoext ])
    out_file = files[b1]
    if (len(out_file)!=1):
        print('WARNING: bin file is repeated or does not exist')
    return (out_file[0])
                  
def run_parspatch(outpath, rhosfile, tausfile,  njk):
    from readjson import read_rhos, read_taus
    from chi2 import minimizeCHI2
    from maxlikelihood import MCMC, percentiles
    import numpy as np

    nwalkers,  nsteps = 100,  1000
    eq = 'All'; moderr = False
    gflag, bflag = True, True
    nsig = 1
    i_guess0 = [ -0.01, 1,- 1] #fiducial values
    data = {}
    data_rhos =  fitsio.read(rhosfile, ext=1)
    data_taus =  fitsio.read(tausfile, ext=1)
    rhosnames =  ['RHO0P','RHO1P','RHO2P','RHO3P','RHO4P','RHO5P']
    varrhosnames =  ['VAR_RHO0','VAR_RHO1','VAR_RHO2','VAR_RHO3','VAR_RHO4','VAR_RHO5']
    tausnames =  ['TAU0P','TAU2P','TAU5P']
    vartausnames =  ['VAR_TAU0','VAR_TAU2','VAR_TAU5']
    for i in range(njk):
        jkrhos = [data_rhos['JKR'] == i]
        jktaus = [data_taus['JKR'] == i]

        a_c = []; a_l = []; a_r = []
        b_c = []; b_l = []; b_r = []
        d_c = []; d_l = []; d_r = []
        for ibin in range(1, 21):
            angrhos = [data_rhos['ANGBIN']< ibin ]
            angtaus = [data_taus['ANGBIN']< ibin ]
            rhosbool = jkrhos&angrhos
            tausbool = jktaus&anftaus
  
            data['rhos'] = [data_rhos[name][rhosbool] for name in rhosnames]
            data['varrhos'] = [data_rhos[name][rhosbool] for name in varrhosnames]
            data['taus'] = [data_taus[name][tausbool] for name in tausnames]
            data['vartaus'] = [data_taus[name][tausbool] for name in vartausnames]

            fit_pars, chisq = minimizeCHI2(data, i_guess0, eq=eq,
                                           gflag=gflag, bflag=bflag,
                                           moderr=moderr)
            samples = MCMC(fit_pars,data, nwalkers, nsteps, eq=eq,
                            gflag=gflag, bflag=bflag, moderr=moderr,
                            plot=False )
            mcmcpars = percentiles(samples, nsig=nsig) 
            a_c.append(mcmcpars[0][0]); a_l.append(mcmcpars[0][1]);a_r.append(mcmcpars[0][2])
            b_c.append(mcmcpars[1][0]); b_l.append(mcmcpars[1][1]);b_r.append(mcmcpars[1][2])
            d_c.append(mcmcpars[2][0]); d_l.append(mcmcpars[2][1]);d_r.append(mcmcpars[2][2])
        write_pars(outpath + 'parspatch' + str(i) + '.json', meanr, a_c, a_l, a_r, b_c, b_l, b_r, d_c, d_l, d_r)
   
def plotlineal(outpath):
    import numpy as np
    colors = ['black', 'green', 'blue', 'red', 'gray', 'pink']
    for patch in range(1, 5):
        meanr, a_c, a_l, a_r, b_c, b_l, b_r, d_c, d_l, d_r = read_pars(outpath+'parspatch'+str(patch)+'.json')
        
        label = "P" + str(patch)
        width = np.diff(meanr)*0.25
        width = np.append(width, width[-1:])
        dx = meanr + (patch - 1)*width
        plt.figure(0)
        plt.scatter(meanr, a_c, color=colors[patch],  label=label, marker='o', s=10)
        plt.errorbar(meanr, a_c, yerr=[-np.array(a_l),np.array(a_r)], color=colors[patch], capsize=2, linestyle='')
        
        plt.figure(1)
        plt.scatter(meanr, b_c, color=colors[patch],  label=label, marker='o', s=10)
        plt.errorbar(meanr, b_c, yerr=[-np.array(b_l),np.array(b_r)], color=colors[patch], capsize=2, linestyle='')
        
        plt.figure(2)
        plt.scatter(meanr, d_c, color=colors[patch],  label=label, marker='o', s=10)
        plt.errorbar(meanr, d_c, yerr=[-np.array(d_l),np.array(d_r)], color=colors[patch], capsize=2, linestyle='')
        
    plt.figure(0)
    plt.legend(loc='best', fontsize=10)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\alpha$', fontsize=24)
    plt.xscale('log')
    plt.xlim( [ 2, 300] )
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
    plt.xlim( [2, 300] )
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
    plt.xlim( [2, 300] )
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
    
    run_parspatch(outpath, args.rhos, args.taus, 4)
    plotlineal(outpath)
    


    
if __name__ == "__main__":
    main()
