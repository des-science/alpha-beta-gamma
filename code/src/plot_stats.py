import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def pretty_rho(meanr, rho, sig,  legend=None, lfontsize=24, color='black', marker='o', ylabel=r'$\rho(\theta)$',title=None,  xlim=None,  ylim=None):
    '''
    plt.plot(meanr, rho, color=color, label=legend, marker=marker)
    plt.plot(meanr, -rho, color=color, ls=':', marker=marker)
    if sig is not None:
        plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color=color, ls='', marker=marker)
        plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color=color, ls='', marker=marker)
        rho0_line = plt.errorbar(-meanr, rho, yerr=sig, color=color, marker=marker)
    '''
    plt.plot(meanr, rho, color=color, label=legend)
    plt.plot(meanr, -rho, color=color, ls=':')
    #rigth quadrants
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color=color, ls='', marker=marker,  capsize=2)
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color=color, ls='', marker=marker,  capsize=2)
    #leftquadrants
    plt.errorbar( -meanr, rho, yerr=sig, color=color,  marker='^',  capsize=2)
    plt.errorbar( -meanr,-rho, yerr=sig, color=color,  marker='^', ls=':', capsize=2)
    plt.legend(loc='best', fontsize=lfontsize)
    if ylim is not None: plt.ylim( ylim )
    if xlim is not None: plt.xlim(xlim)
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(ylabel, fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    if(title): plt.title(title)
    plt.tight_layout()

def pretty_rho0(meanr, rho, sig, title= None,  xlim=None, ylim=None):
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color='blue', ls='', marker='o')
    rho0_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    rho0_line = plt.errorbar(-meanr,- rho, yerr=sig, color='blue', marker='o', ls=':')
    plt.legend([rho0_line],[r'$\rho_0(\theta)$'],loc='upper right', fontsize=24)
    
    plt.tick_params(axis='both', which='major', labelsize=24)
    if ylim is not None: plt.ylim( ylim )
    if xlim is not None: plt.xlim(xlim)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\rho(\theta)$', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    if(title): plt.title(title)
    plt.tight_layout()

def pretty_rho1(meanr, rho, sig, rho3=None, sig3=None, rho4=None, sig4=None, title= None,xlim=None, ylim=None):
    import matplotlib.patches as mp
    
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color='blue', ls='', marker='o')
    rho1_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    #rho1_line = plt.errorbar(-meanr,-rho, yerr=sig, color='blue', marker='o', ls=':')
    if rho3 is not None:
        plt.plot(meanr*1.03, rho3, color='green')
        plt.plot(meanr*1.03, -rho3, color='green', ls=':')
        plt.errorbar(meanr[rho3>0]*1.03, rho3[rho3>0], yerr=sig3[rho3>0], color='green', ls='', marker='s')
        plt.errorbar(meanr[rho3<0]*1.03, -rho3[rho3<0], yerr=sig3[rho3<0], color='green', ls='', marker='s')
        rho3_line = plt.errorbar(-meanr, rho3, yerr=sig3, color='green', marker='s')
        #rho3_line = plt.errorbar(-meanr,-rho3, yerr=sig3, color='green', marker='s', ls=':')
    if rho4 is not None:
        plt.plot(meanr*1.06, rho4, color='red')
        plt.plot(meanr*1.06, -rho4, color='red', ls=':')
        plt.errorbar(meanr[rho4>0]*1.06, rho4[rho4>0], yerr=sig4[rho4>0], color='red', ls='', marker='^')
        plt.errorbar(meanr[rho4<0]*1.06, -rho4[rho4<0], yerr=sig4[rho4<0], color='red', ls='', marker='^')
        rho4_line = plt.errorbar(-meanr, rho4, yerr=sig4, color='red', marker='^')
        #rho4_line = plt.errorbar(-meanr,-rho4, yerr=sig4, color='red', marker='^', ls=':')

    plt.legend([rho1_line, rho3_line, rho4_line],
               [r'$\rho_1(\theta)$', r'$\rho_3(\theta)$',
                r'$\rho_4(\theta)$'], loc='upper right',
               fontsize=24)

    plt.tick_params(axis='both', which='major', labelsize=24)
    if ylim is not None: plt.ylim( ylim )
    if xlim is not None: plt.xlim(xlim)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\rho(\theta)$', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    if(title): plt.title(title)
    plt.tight_layout()

def pretty_rho2(meanr, rho, sig,  rho5=None, sig5=None, title= None, xlim=None, ylim=None ):
    import matplotlib.patches as mp
    # The requirements on rho2 are less stringent.  They are larger by a factor 1/alpha.
    # Let's use alpha = 0.03.
    alpha = 0.03

    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color='blue', ls='', marker='o')
    rho2_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    #rho2_line = plt.errorbar(-meanr,-rho, yerr=sig, color='blue', marker='o', ls=':')
    if rho5 is not None:
        plt.plot(meanr*1.03, rho5, color='green')
        plt.plot(meanr*1.03, -rho5, color='green', ls=':')
        plt.errorbar(meanr[rho5>0]*1.03, rho5[rho5>0], yerr=sig5[rho5>0], color='green', ls='', marker='s')
        plt.errorbar(meanr[rho5<0]*1.03, -rho5[rho5<0], yerr=sig5[rho5<0], color='green', ls='', marker='s')
        rho5_line = plt.errorbar(-meanr, rho5, yerr=sig5, color='green', marker='s')
        #rho5_line = plt.errorbar(-meanr,-rho5, yerr=sig5, color='green', marker='s', ls=':')
        
    if rho5 is not None :
        plt.legend([rho2_line, rho5_line ],
                   [r'$\rho_2(\theta)$', r'$\rho_5(\theta)$'],
                   loc='upper right', fontsize=24)

    plt.tick_params(axis='both', which='major', labelsize=24)
    if ylim is not None: plt.ylim( ylim )
    if xlim is not None: plt.xlim(xlim)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\rho(\theta)$', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    if(title): plt.title(title)
    plt.tight_layout()

def pretty_tau(meanr, rho, sig, title= None, xlim=None,  ylim=None):
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0], color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0], color='blue', ls='', marker='o')
    rho0_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    #rho0_line = plt.errorbar(-meanr,-rho, yerr=sig, color='blue', marker='o', ls=':')
    plt.legend([rho0_line],[title],loc='upper right', fontsize=24)
    #plt.ylim( [1.e-9, 5.e-6] )
    plt.tick_params(axis='both', which='major', labelsize=24)
    if ylim is not None: plt.ylim( ylim )
    if xlim is not None: plt.xlim(xlim)
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(title, fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    if(title): plt.title(title)
    plt.tight_layout()

def plotcorrmat(cov):
    import numpy as np
    cov = np.mat(cov)
    D = np.diag(np.sqrt(np.diag(cov)))
    d = np.linalg.inv(D)
    corr = d*cov*d
    cov_vmin=np.min(corr)
    plt.imshow(corr,cmap='viridis'+'_r', interpolation='nearest',
               aspect='auto', origin='lower', vmin=cov_vmin, vmax=1.)
    plt.colorbar()
    
def plotallrhos(filename, outpath, title= None, xlim=None, ylims=None):
    from readjson import read_rhos
    ylim0, ylim1, ylim2 = [None, None, None]
    if ylims is not None: ylim0, ylim1, ylim2 = ylims
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 = read_rhos(filename)
    plt.clf()
    pretty_rho0(meanr, rho0p, sig_rho0, title=title, xlim=xlim, ylim=ylim0)
    print("Printing file: ", outpath +'rho0_all_rsrs.png')
    plt.savefig(outpath +'rho0_all_rsrs.png')
    plt.clf()
    pretty_rho1(meanr, rho1p, sig_rho1,  rho3p, sig_rho3, rho4p, sig_rho4, title=title, xlim=xlim,  ylim=ylim1)
    print("Printing file: ", outpath +'rho1_all_rsrs.png')
    plt.savefig(outpath +'rho1_all_rsrs.png')
    plt.clf()
    pretty_rho2(meanr, rho2p, sig_rho2,  rho5p, sig_rho5, title=title, xlim=xlim,  ylim=ylim2)
    print("Printing file: ", outpath +'rho2_all_rsrs.png')
    plt.savefig(outpath +'rho2_all_rsrs.png')
    
def plotallrhosfits(filenames, outpath, title= None, xlim=None, ylims=None):
    import numpy as np
    from readfits import read_corr
    ylim0, ylim1, ylim2 = [None, None, None]
    if ylims is not None: ylim0, ylim1, ylim2 = ylims
    meanr, rho0, cov_rho0 = read_corr(filenames[0])
    meanr, rho1, cov_rho1 = read_corr(filenames[1])
    meanr, rho2, cov_rho2 = read_corr(filenames[2])
    meanr, rho3, cov_rho3 = read_corr(filenames[3])
    meanr, rho4, cov_rho4 = read_corr(filenames[4])
    meanr, rho5, cov_rho5 = read_corr(filenames[5])
    sig_rho0 =  np.sqrt(np.diag(cov_rho0))
    sig_rho1 =  np.sqrt(np.diag(cov_rho1))
    sig_rho2 =  np.sqrt(np.diag(cov_rho2))
    sig_rho3 =  np.sqrt(np.diag(cov_rho3))
    sig_rho4 =  np.sqrt(np.diag(cov_rho4))
    sig_rho5 =  np.sqrt(np.diag(cov_rho5))
    plt.clf()
    pretty_rho0(meanr, rho0, sig_rho0, title=title, xlim=xlim, ylim=ylim0)
    print("Printing file: ", outpath +'rho0_all_rsrs.png')
    plt.savefig(outpath +'rho0_all_rsrs.png')
    plt.clf()
    pretty_rho1(meanr, rho1, sig_rho1, rho3, sig_rho3, rho4, sig_rho4, title=title, xlim=xlim, ylim=ylim1)
    print("Printing file: ", outpath +'rho1_all_rsrs.png')
    plt.savefig(outpath +'rho1_all_rsrs.png')
    plt.clf()
    pretty_rho2(meanr, rho2, sig_rho2,  rho5, sig_rho5, title=title, xlim=xlim, ylim=ylim2)
    print("Printing file: ", outpath +'rho2_all_rsrs.png')
    plt.savefig(outpath +'rho2_all_rsrs.png')
    
def plotalltausfits(filenames, outpath, title= None,xlim=None,  ylims=None):
    import numpy as np
    from readfits import read_corr
    ylim0, ylim1, ylim2 = [None, None, None]
    if ylims is not None: ylim0, ylim1, ylim2 = ylims
    meanr, tau0, cov_tau0 = read_corr(filenames[0])
    meanr, tau2, cov_tau2 = read_corr(filenames[1])
    meanr, tau5, cov_tau5 = read_corr(filenames[2])
    sig_tau0 =  np.sqrt(np.diag(cov_tau0))
    sig_tau2 =  np.sqrt(np.diag(cov_tau2))
    sig_tau5 =  np.sqrt(np.diag(cov_tau5))
    plt.clf()
    pretty_tau(meanr2, tau0, sig_tau0, r'$\tau_{0}(\theta)$',title=title,  xlim=xlim, ylim=ylim0)
    print("Printing file: ", outpath +'tau0_all_rsgal.png')
    plt.savefig(outpath +'tau0_all_rsgal.png')
    plt.clf()
    pretty_tau(meanr2, tau2, sig_tau2, r'$\tau_{2}(\theta)$',title=title,  xlim=xlim, ylim=ylim1)
    print("Printing file: ", outpath +'tau2_all_rsgal.png')
    plt.savefig(outpath +'tau2_all_rsgal.png')
    plt.clf()
    pretty_tau(meanr2, tau5, sig_tau5, r'$\tau_{5}(\theta)$',title=title,  xlim=xlim, ylim=ylim2)
    print("Printing file: ", outpath +'tau5_all_rsgal.png')
    plt.savefig(outpath +'tau5_all_rsgal.png')
def plotallrhoscorrmatfits(filenames, outpath):
    import fitsio
    names = ['rho0_covmat.png', 'rho1_covmat.png', 'rho2_covmat.png', 'rho3_covmat.png', 'rho4_covmat.png', 'rho5_covmat.png']
    titles =  [r'$\rho_{0}(\theta)$',r'$\rho_{1}(\theta)$', r'$\rho_{2}(\theta)$', r'$\rho_{3}(\theta)$', r'$\rho_{4}(\theta)$', r'$\rho_{5}(\theta)$' ]
    for i,f in enumerate(filenames):
        covmat = fitsio.read(f, ext=1)
        plt.clf()
        plotcorrmat(covmat)
        plt.title(titles[i])
        plt.savefig(outpath + names[i], dpi=500)
        print(outpath +names[i], 'Printed!')
def plotalltauscorrmatfits(filenames, outpath):
    import fitsio
    names = ['tau0_covmat.png', 'tau2_covmat.png', 'tau_covmat.png']
    titles =  [r'$\rho_{0}(\theta)$', r'$\rho_{2}(\theta)$', r'$\rho_{5}(\theta)$' ]
    for i,f in enumerate(filenames):
        covmat = fitsio.read(f, ext=1)
        plt.clf()
        plotcorrmat(covmat)
        plt.title(titles[i])
        plt.savefig(outpath + names[i], dpi=500)
        print(outpath +names[i], 'Printed!')
def plotalltaus(filename, outpath, xlim=None, ylims=None):
    from readjson import read_taus
    ylim0, ylim1, ylim2 = [None, None, None]
    if ylims is not None: ylim0, ylim1, ylim2 = ylims
    meanr2, tau0p, tau2p, tau5p, sig_tau0, sig_tau2, sig_tau5 =  read_taus(filename)
    plt.clf()
    pretty_tau(meanr2, tau0p, sig_tau0, r'$\tau_{0}(\theta)$', xlim=xlim, ylim=ylim0)
    print("Printing file: ", outpath +'tau0_all_rsgal.png')
    plt.savefig(outpath +'tau0_all_rsgal.png')
    plt.clf()
    pretty_tau(meanr2, tau2p, sig_tau2, r'$\tau_{2}(\theta)$', xlim=xlim, ylim=ylim1)
    print("Printing file: ", outpath +'tau2_all_rsgal.png')
    plt.savefig(outpath +'tau2_all_rsgal.png')
    plt.clf()
    pretty_tau(meanr2, tau5p, sig_tau5, r'$\tau_{5}(\theta)$', xlim=xlim, ylim=ylim2)
    print("Printing file: ", outpath +'tau5_all_rsgal.png')
    plt.savefig(outpath +'tau5_all_rsgal.png')

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
  
def plot_tomograpically(ax, x=None, ypath=None,  xerr=None, ypatherr=None,xlabel='', ylabel='', nbins=None, nbins1=None, nbins2=None, symm=False, color='blue'):
    # This function is defined to work when the output are in the
    # usual form of comosis, i.e one common x file, a possible x error
    # file. And multiples ipath(yfiles), using the format name_i_j.txt
    import itertools
    import os
    if(nbins):
        nbins1=nbins;nbins2=nbins
        symm=True
    elif(nbins1 and nbins2):
        if(symm==True):
            print("Warning, Symmetry requested by differenve nbins")
        symm=False
    else:
        print('Error, defining symultaneosly nbins and nbins1 or nbins2')
        raise
        
    a=[i for i in range(1,nbins1+1)]
    b=[j for j in range(1,nbins2+1)]
    bin_pairs=[]
    for p in itertools.product(a, b):
        bin_pairs.append(p)
    mtype = '.'
    
    for i,j in bin_pairs:
        file = findbinfile(ypath, i, j , symm=symm)
        yname= os.path.join(ypath, file);y = np.loadtxt(yname)
        if (ypatherr):
            fileerr = findbinfile_sym(ypatherr, i, j )
            yerrname= os.path.join(ypatherr, fileerr); yerr = np.loadtxt(yerrname)
        else:
            yerror = np.zeros(len(y))
        if (xerr):
            fileerr = findbinfile_sym(xpatherr, i, j )
            xerrname= os.path.join(ypatherr, fileerr); xerr = np.loadtxt(yerrname)
        else:
            yerror = np.zeros(len())
        ax[j-1][i-1].errorbar(x, abs(y),xerr=xerror,yerr=yerror, fmt=mtype ,capsize=0,markersize=1, color=color, mec=color, elinewidth=1.)
        ax[j-1][i-1].text(0.85, 0.85, "{},{}".format(i, j), horizontalalignment='center', verticalalignment='center', transform=ax[j-1][i-1].transAxes, fontsize=12)
        ax[j-1][i-1].set_xscale('log', nonposx='clip')
        ax[j-1][i-1].set_yscale('log', nonposy='clip')
        if (j == nbins2):
            ax[j-1][i-1].set_xlabel(xlabel)
        if (i == 1):
            ax[j-1][i-1].set_ylabel(ylabel)
        if (i>j and symm):
            ax[j-1][i-1].set_visible(False)

def plot_tomograpically_bin(ax, i, j,  x, y, xerr=None, yerr=None,xlabel='', ylabel='', nbins=None, nbins1=None, nbins2=None, color='blue',label=''):
    import numpy as np
    if(nbins):
        nbins1=nbins;nbins2=nbins
    elif(nbins1 and nbins2):
        print("Entering in not symmetric mode")
    else:
        print('Error, defining symultaneosly nbins and nbins1 or nbins2')
        raise
        
    mtype = '.'
    if xerr is None: xerr =  np.zeros(len(x))
    if yerr is None: yerr =  np.zeros(len(y))
    ax[j-1][i-1].errorbar(x, abs(y),xerr=xerr,yerr=yerr, fmt=mtype ,capsize=0,markersize=1, color=color, mec=color, elinewidth=0.5, label=label)
    ax[j-1][i-1].text(0.85, 0.85, "{},{}".format(i, j), horizontalalignment='center', verticalalignment='center', transform=ax[j-1][i-1].transAxes, fontsize=12)
    ax[j-1][i-1].set_xscale('log', nonposx='clip')
    ax[j-1][i-1].set_yscale('log', nonposy='clip')
    ax[j-1][i-1].legend(loc='best')
    if (j == nbins2):
        ax[j-1][i-1].set_xlabel(xlabel)
    if (i == 1):
        ax[j-1][i-1].set_ylabel(ylabel)
