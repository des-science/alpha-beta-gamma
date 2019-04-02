#Ploting all tau
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Plot to compare xipobs with possible due to psf modelling ')
    
    parser.add_argument('--tausbins',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/tomotaus/',
                        help='Json file with the xip_obs correlation function')
    parser.add_argument('--srcpath',
                         default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src',
                         help='Path to src/ folder with plenty functions required')
    parser.add_argument('--outpath',
                        default='/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/plots/',
                        help='location of the output of the files, by default plots is created')
    
    
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
            print('WARNING: bin file is repeated, probably you 2pt is not simmetric')
        return (out_file[0])
    else:
        out_file = files[b1]
        if (len(out_file)!=1):
            print('WARNING: bin file is repeated or does not exist')
        return (out_file[0])
  
def main():
    import sys
    args = parse_args()
    sys.path.insert(0, args.srcpath)
    import numpy as np
    from readjson import read_taus
    from plot_stats import pretty_rho, plot_tomograpically_bin

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    section_name= 'galaxy_cl_nbb'
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\tau$'
    nbins=4
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    bins = [[1, 1], [2, 2], [3, 3], [4, 4]]
    lines = []; cols = ['red', 'blue', 'green', 'gray']
    for i, j in bins:
        fileaux = findbinfile(args.tausbins, i, j)
        meanr, tau0p, tau2p, tau5p, sig_tau0p, sig_tau2p, sig_tau5p = read_taus(args.tausbins + fileaux)
        plot_tomograpically_bin(ax, i, j, meanr, tau0p,
                                yerr=sig_tau0p, xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue')#, label=r'$\tau_{0}$')
        plot_tomograpically_bin(ax, i, j, meanr, tau2p,
                                yerr=sig_tau2p, xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='red')#, label=r'$\tau_{2}$')
        plot_tomograpically_bin(ax, i, j, meanr, tau5p,
                                yerr=sig_tau5p, xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='green')#, label=r'$\tau_{05}$')
        lines.append(ax[1][1].lines)

    fig.legend(lines,labels=[r'$\tau_{0}$', r'$\tau_{2}$', r'$\tau_{5}$'  ], bbox_to_anchor=(1.1, 1.08))
    plt.tight_layout()
    fname = 'all_taus_tomo.png'
    print('Printing:',  outpath +fname)
    plt.savefig(outpath +fname,  dpi=300)
    plt.close(fig)

    #tau0
    plt.clf()
    for i, j in bins:
        fileaux = findbinfile(args.tausbins, i, j)
        meanr, tau0p, tau2p, tau5p, sig_tau0p, sig_tau2p, sig_tau5p = read_taus(args.tausbins + fileaux)
        pretty_rho(meanr, tau0p, sig_tau0p, legend="{},{}".format(i, j), lfontsize=15,  color=cols[i-1],  ylabel=r'$\tau_{0}(\theta)$', ylim=False)
    fname = 'taus0_tomo.png'
    print('Printing:',  outpath +fname)
    plt.savefig(outpath +fname,  dpi=300)
    #tau2
    plt.clf()
    for i, j in bins:
        fileaux = findbinfile(args.tausbins, i, j)
        meanr, tau0p, tau2p, tau5p, sig_tau0p, sig_tau2p, sig_tau5p = read_taus(args.tausbins + fileaux)
        pretty_rho(meanr, tau2p, sig_tau2p, legend="{},{}".format(i, j),lfontsize=15, color=cols[i-1],  ylabel=r'$\tau_{2}(\theta)$', ylim=False)
    fname = 'taus2_tomo.png'
    print('Printing:',  outpath +fname)
    plt.savefig(outpath +fname,  dpi=300)
    #tau5
    plt.clf()
    for i, j in bins:
        fileaux = findbinfile(args.tausbins, i, j)
        meanr, tau0p, tau2p, tau5p, sig_tau0p, sig_tau2p, sig_tau5p = read_taus(args.tausbins + fileaux)
        pretty_rho(meanr, tau5p, sig_tau5p, legend="{},{}".format(i, j), lfontsize=15,color=cols[i-1],  ylabel=r'$\tau_{5}(\theta)$', ylim=False)
    fname = 'taus5_tomo.png'
    print('Printing:',  outpath +fname)
    plt.savefig(outpath +fname,  dpi=300)




    
    

if __name__ == "__main__":
    main()
    
