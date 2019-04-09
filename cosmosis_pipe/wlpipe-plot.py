import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#plt.style.use('/data/git_repositories/alpha-beta-gamma/code/SVA1StyleSheet.mplstyle')
import matplotlib.colors as colors


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Sofware to plot all quantities after running WL-pipeline')
    
    parser.add_argument('--fiducial',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/cosmosis_pipe/2pt_sim_1110_baseline_Y3cov.fits',
                        help='fit file with fiducial data vectors, covariance matrix and so on.')
    parser.add_argument('--contaminant',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/cosmosis_pipe/abg_dxip_tomo.fits',
                        help='fit file with contamination data vector, covariance matrix')
    parser.add_argument('--contaminated',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/cosmosis_pipe/2pt_sim_1110_baseline_Y3cov_contaminated.fits', 
                        help='fit file with contamination data vector, covariance matrix')
    parser.add_argument('--out',
                        default='/home/dfa/sobreira/alsina/alpha-beta-gamma/cosmosis_pipe/',
                        help='path where output will be send') 
    
    args = parser.parse_args()
    return args

def plotfiducial(fitfile, out):
    import fitsio
    import itertools
    import numpy as np
    #PLOTING FIDUCIAL INFO
    fiducialfit = fitfile
    covmatrixfit=fitsio.read(fiducialfit,ext=1)
    xipfit=fitsio.read(fiducialfit,ext=2)
    ximfit=fitsio.read(fiducialfit,ext=3)
    gammatfit=fitsio.read(fiducialfit,ext=4)
    wthetafit=fitsio.read(fiducialfit,ext=5)
    nz_sourcefit=fitsio.read(fiducialfit,ext=6)
    nz_lensfit=fitsio.read(fiducialfit,ext=7)
    lengths = [len(xipfit),len(ximfit),len(gammatfit),len(wthetafit)]

    ##Covariance Matrix.
    corr=corrmatrix(covmatrixfit)
    cov_vmin=np.min(corr)
    plt.imshow(corr,cmap='viridis'+'_r', interpolation='nearest',
               aspect='auto', origin='lower', vmin=cov_vmin, vmax=1.)
    plt.colorbar()
    plt.title(r'$\xi_{+}(\theta) \mid \xi_{-}(\theta) \mid \gamma_{t}(\theta) \mid \omega(\theta)$')
    pos_lines = [0]
    for i in range(len(lengths)):
        pos_lines.append(pos_lines[i] + lengths[i])
    pos_lines = pos_lines[1:-1]
    for line in pos_lines:
        plt.axvline(x=line, c='k', lw=1, ls='-')
        plt.axhline(y=line, c='k', lw=1, ls='-')
    plt.tight_layout()
    filename = out + 'CovariancematrixFiducial.png'
    plt.savefig(filename, dpi=500)
    print(filename, 'Printed!')

    #xip
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{+}(\theta)$'
    nbins=4
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin = (xipfit['BIN1']==i)&(xipfit['BIN2']==j)
        theta=xipfit['ANG'][bin]
        xip=xipfit['VALUE'][bin]
        yerr=get_error(covmatrixfit, lengths, 'xip')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xip, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue')
    filename = out + 'xip_fiducial.png'
    plt.savefig(filename, dpi=300)
    print(filename, 'Printed!')
    
    ##xim
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{-}(\theta)$'
    nbins=4
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin = (ximfit['BIN1']==i)&(ximfit['BIN2']==j)
        theta=ximfit['ANG'][bin]
        xim=ximfit['VALUE'][bin]
        yerr=get_error(covmatrixfit, lengths, 'xim')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xim, yerr=yerr, xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue')
    filename = out + 'xim_fiducial.png'
    plt.savefig(filename, dpi=300)
    print(filename, 'Printed!')
    
    ##wtheta
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\omega(\theta)$'
    nbins=5
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin = (wthetafit['BIN1']==i)&(wthetafit['BIN2']==j)
        theta=wthetafit['ANG'][bin]
        wtheta=wthetafit['VALUE'][bin]
        yerr=get_error(covmatrixfit, lengths, 'wtheta')[bin]
        plot_tomograpically_bin(ax, i, j, theta, wtheta, yerr=yerr, xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='darkblue')
    plt.tight_layout()
    filename = out + 'wtheta_fiducial.png'
    plt.savefig(filename, dpi=300)
    print(filename, 'Printed!')
def plotcontaminant(fitfile, out):
    import fitsio
    import itertools
    import numpy as np
    contaminantfit = fitfile
    covmatrixfit=fitsio.read(contaminantfit,ext=1)
    xipfit=fitsio.read(contaminantfit,ext=2)
    lengths = [len(xipfit)]
    
    ##Covariance Matrix.
    corr=corrmatrix(covmatrixfit)
    cov_vmin=np.min(corr)
    plt.imshow(corr,cmap='viridis'+'_r', interpolation='nearest',
               aspect='auto', origin='lower', vmin=cov_vmin, vmax=1.)
    plt.colorbar()
    plt.title(r'$\xi_{+}(\theta)$')
    pos_lines = [0]
    for i in range(len(lengths)):
        pos_lines.append(pos_lines[i] + lengths[i])
    pos_lines = pos_lines[1:-1]
    for line in pos_lines:
        plt.axvline(x=line, c='k', lw=1, ls='-')
        plt.axhline(y=line, c='k', lw=1, ls='-')
    plt.tight_layout()
    filename = out + 'Covariancematrix_contaminant.png'
    plt.savefig(filename, dpi=500)
    print(filename, 'Printed!')

    #xip
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\delta \xi_{+}(\theta)$'
    nbins=4
    fig, ax = plt.subplots(1, nbins, figsize=(1.6*nbins, 1.6), sharey=True, sharex=True)
    bin_pairs=[[1, 1], [2, 2], [3, 3], [4, 4]]
    for i,j in bin_pairs:
        bin = (xipfit['BIN1']==i)&(xipfit['BIN2']==j)
        theta=xipfit['ANG'][bin]
        xip=xipfit['VALUE'][bin]
        yerr=get_error(covmatrixfit, lengths, 'xip')[bin]
        ax[i-1].errorbar(theta,
                              abs(xip),xerr=np.zeros(len(theta)),yerr=yerr,
                              fmt='o' ,capsize=1,markersize=2,
                              color='red', mec='red', elinewidth=0.5,
                              label='')
        ax[i-1].text(0.85, 0.85, "{},{}".format(i, j), horizontalalignment='center', verticalalignment='center', transform=ax[i-1].transAxes, fontsize=12)
        ax[i-1].set_xscale('log', nonposx='clip')
        ax[i-1].set_yscale('log', nonposy='clip')

    filename = out + 'xip_contaminant.png'
    plt.savefig(filename, dpi=300)
    print(filename, 'Printed!')
    
def plotcontaminated(fitfile, out):
    import fitsio
    import itertools
    import numpy as np
    #PLOTING FIDUCIAL INFO
    fiducialfit = fitfile
    covmatrixfit=fitsio.read(fiducialfit,ext=1)
    xipfit=fitsio.read(fiducialfit,ext=2)
    ximfit=fitsio.read(fiducialfit,ext=3)
    gammatfit=fitsio.read(fiducialfit,ext=4)
    wthetafit=fitsio.read(fiducialfit,ext=5)
    nz_sourcefit=fitsio.read(fiducialfit,ext=6)
    nz_lensfit=fitsio.read(fiducialfit,ext=7)
    lengths = [len(xipfit),len(ximfit),len(gammatfit),len(wthetafit)]
    
    ##Covariance Matrix.
    corr=corrmatrix(covmatrixfit)
    cov_vmin=np.min(corr)
    plt.imshow(corr,cmap='viridis'+'_r', interpolation='nearest',
               aspect='auto', origin='lower', vmin=cov_vmin, vmax=1.)
    plt.colorbar()
    plt.title(r'$\xi_{+}(\theta) \mid \xi_{-}(\theta) \mid \gamma_{t}(\theta) \mid \omega(\theta)$')
    pos_lines = [0]
    for i in range(len(lengths)):
        pos_lines.append(pos_lines[i] + lengths[i])
    pos_lines = pos_lines[1:-1]
    for line in pos_lines:
        plt.axvline(x=line, c='k', lw=1, ls='-')
        plt.axhline(y=line, c='k', lw=1, ls='-')
    plt.tight_layout()
    filename = out + 'Covariancematrixcontaminated.png'
    plt.savefig(filename, dpi=500)
    print(filename, 'Printed!')

    #xip
    xlabel=r'$\theta$ [arcmin]'; ylabel=r'$\xi_{+}(\theta)$'
    nbins=4
    fig, ax = plt.subplots(nbins, nbins, figsize=(1.6*nbins, 1.6*nbins), sharey=True, sharex=True)
    a=[i for i in range(1,nbins+1)]
    b=[j for j in range(1,nbins+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    for i,j in bin_pairs:
        bin = (xipfit['BIN1']==i)&(xipfit['BIN2']==j)
        theta=xipfit['ANG'][bin]
        xip=xipfit['VALUE'][bin]
        yerr=get_error(covmatrixfit, lengths, 'xip')[bin]
        plot_tomograpically_bin(ax, i, j, theta, xip, yerr=yerr,
                                xlabel=xlabel, ylabel=ylabel, nbins=4,
                                color='blue')
    filename = out + 'xip_contaminated.png'
    plt.savefig(filename, dpi=300)
    print(filename, 'Printed!')
def checkcontamination():
def plot_tomograpically_bin(ax, i, j, x, y, xerr=None,
                            yerr=None,xlabel='', ylabel='',
                            nbins=None, nbins1=None, nbins2=None,
                            color='blue',label=None):
    import numpy as np
    if(nbins):
        nbins1=nbins;nbins2=nbins
    elif(nbins1 and nbins2):
        print("Entering in not symmetric mode")
    else:
        print('Error, defining symultaneosly nbins and nbins1 or nbins2')
        raise
        
    mtype = 'o'
    if xerr is None: xerr =  np.zeros(len(x))
    if yerr is None: yerr =  np.zeros(len(y))
    
    if y.size !=0 :
        ax[j-1][i-1].errorbar(x, abs(y),xerr=xerr,yerr=yerr, fmt=mtype ,capsize=1,markersize=2, color=color, mec=color, elinewidth=0.5, label=label)
        ax[j-1][i-1].text(0.85, 0.85, "{},{}".format(i, j), horizontalalignment='center', verticalalignment='center', transform=ax[j-1][i-1].transAxes, fontsize=12)
        ax[j-1][i-1].set_xscale('log', nonposx='clip')
        ax[j-1][i-1].set_yscale('log', nonposy='clip')
        if(label):
            ax[j-1][i-1].legend(loc='best')
        if (j == nbins2):
            ax[j-1][i-1].set_xlabel(xlabel)
        if (i == 1):
            ax[j-1][i-1].set_ylabel(ylabel)
    else:
        ax[j-1][i-1].set_visible(False)
def get_error(covmatrix, lengths, name):
    import numpy as np
    if name is not None:
        if (name=='xip'):
            start = 0
            end =start + lengths[0]
        elif (name=='xim'):
            start = lengths[0]
            end =start + lengths[1]
        elif (name=='gammat'):
            start = lengths[0] + lengths[1]
            end =start + lengths[2]
        elif (name=='wtheta'):
            start = lengths[0] + lengths[1]+ lengths[2]
            end =start + lengths[3]
        return np.diagonal(covmatrix)[start:end]**0.5
    else:
        print("Correlation function no defined")
        return None

##Ploting covariance matrix might be not than ilustrative as
##correlations matrix
def corrmatrix(cov):
    import numpy as np
    cov = np.mat(cov)
    D = np.diag(np.sqrt(np.diag(cov)))
    d = np.linalg.inv(D)
    corr = d*cov*d
    return corr

def main():
    
    import numpy as np
    args = parse_args()
    out = os.path.expanduser(args.out + 'plots/')
    try:
        if not os.path.isdir(out):
            os.makedirs(out)
    except OSError as e:
        print "Ignore OSError from makedirs(work):"
        print e
        pass
    
    #plotfiducial(args.fiducial,  out)
    #plotcontaminant(args.contaminant, out)
    plotcontaminated(args.contaminated, out)
    


if __name__ == "__main__":
    main()
