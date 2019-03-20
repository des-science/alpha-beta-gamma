import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('/data/git_repositories/alpha-beta-gamma/code/SVA1StyleSheet.mplstyle')
from matplotlib.pyplot import figure
#figure(num=None, figsize=(10, 8), dpi=150, facecolor='w', edgecolor='k')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Sofware to plot all quantities after running WL-pipeline')
    
    parser.add_argument('--work',
                        default='/data/git_repositories/alpha-beta-gamma/cosmosis_pipe/test_output_d_l',
                        help='Master path with all the folder generated by the pipeline')
    parser.add_argument('--out',
                        default='/data/git_repositories/alpha-beta-gamma/cosmosis_pipe/',
                        help='Path for the outputs of this code')  
    
    args = parser.parse_args()
    return args

def plotTGraph2D(data, field, xtitle, ytitle, outfile_name):
    from ROOT import TCanvas, TGraph,  TH2F,  TH2,  TH1,  TH1F
    from ROOT import gROOT, gSystem,  Double,  TAxis,  gStyle
    import numpy as np

    ignore =  (data[field] ==-999.) |  (data[field] == None) | ( np.isnan(data[field])) 
    data = data[~ignore]
    
    gStyle.SetOptStat(0)
    c1 =  TCanvas('c1', '', 1000,1000)
    c1.SetBottomMargin( 0.15 )
    c1.SetTopMargin( 0.05 )
    c1.SetLeftMargin( 0.15 )
    c1.SetRightMargin( 0.15 )
    c1.Divide(1, 2)

    data = data[(data['rho2p']>0)]
   
    g = TGraph(len(data[field]), data[field],  data['rho2p'] )
   
    c1.cd().SetLogy()
    g.Draw('AC*')
    g.SetTitleOffset(  0.6 , "y"); 
    g.SetTitleOffset(  0.6 , "x"); 
    g.GetYaxis().CenterTitle()
    g.GetYaxis().SetTitle(ytitle)
    g.GetXaxis().CenterTitle()
    g.GetXaxis().SetTitle(xtitle)
    g.SetTitleSize(0.065, "x"); 
    g.SetTitleSize(0.065, "y");  
    c1.Print(outfile_name)
def camb(cambpath,  out):
    import numpy as np
    linewidth=4
    fontsize = 34

    # CMB
    cmbpath =  os.path.join(cambpath, 'cmb_cl/')
    ellname = os.path.join(cmbpath, 'ell.txt'); ell = np.loadtxt(ellname)
    bbname = os.path.join(cmbpath,  'bb.txt'); bb = np.loadtxt(bbname)
    eename = os.path.join(cmbpath,  'ee.txt'); ee = np.loadtxt(eename)
    tename = os.path.join(cmbpath,  'te.txt'); te = np.loadtxt(tename)
    ttname = os.path.join(cmbpath,  'tt.txt'); tt = np.loadtxt(ttname)

    plt.clf()
    plt.xticks(fontsize=fontsize); plt.yticks(fontsize=fontsize)
    plt.plot(ell, bb,  linewidth=linewidth)
    plt.xlabel(r'l', fontsize=34)
    plt.ylabel(r'l(l+1)C$_{l}$/2$\pi$BB/$\mu$K$^{2}$', fontsize=24)
    plt.xlim([0, 400]);  #plt.ylim([-260, 260])
    name = out + 'bb.png'
    plt.tight_layout()
    plt.savefig(name)
    print(name, 'PRINTED!')

    plt.clf()
    plt.xticks(fontsize=fontsize); plt.yticks(fontsize=fontsize)
    plt.plot(ell, ee, linewidth=linewidth)
    plt.xlabel(r'l', fontsize=fontsize)
    plt.ylabel(r'l(l+1)C$_{l}$/2$\pi$EE/$\mu$K$^{2}$', fontsize=24)
    #plt.xlim([0, 400]);  #plt.ylim([-260, 260])
    name = out + 'ee.png'
    plt.tight_layout()
    plt.savefig(name)
    print(name, 'PRINTED!')

    plt.clf()
    plt.xticks(fontsize=fontsize); plt.yticks(fontsize=fontsize)
    plt.plot(ell, te, linewidth=linewidth)
    plt.xlabel(r'l', fontsize=fontsize)
    plt.ylabel(r'l(l+1)C$_{l}$/2$\pi$TE/$\mu$K$^{2}$', fontsize=24)
    #plt.xlim([0, 400]);  #plt.ylim([-260, 260])
    name = out + 'te.png'
    plt.tight_layout()
    plt.savefig(name)
    print(name, 'PRINTED!')

    plt.clf()
    plt.xticks(fontsize=fontsize); plt.yticks(fontsize=fontsize)
    plt.plot(ell, tt, linewidth=linewidth)
    plt.xlabel(r'l', fontsize=fontsize)
    plt.ylabel(r'l(l+1)C$_{l}$/2$\pi$TT/$\mu$K$^{2}$', fontsize=24)
    #plt.xlim([0, 400]);  #plt.ylim([-260, 260])
    name = out + 'tt.png'
    plt.tight_layout()
    plt.savefig(name)
    print(name, 'PRINTED!')


    #TRANSFER FUNCTION
    linear_cdm_transfer_path =  os.path.join(cambpath, 'linear_cdm_transfer/')
    zname = os.path.join(linear_cdm_transfer_path, 'z.txt'); z = np.loadtxt(zname)
    khname = os.path.join(linear_cdm_transfer_path, 'k_h.txt'); kh = np.loadtxt(khname)
    deltaname = os.path.join(linear_cdm_transfer_path, 'delta_cdm.txt'); delta = np.loadtxt(deltaname)
    
    X, Y = np.meshgrid(kh,z); Z = delta
    plt.clf()
    fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    ax = Axes3D(fig)
    ax.plot_surface(X, Y, Z, linewidth=0, rstride=4, cstride=4)
    #ax.contour(X, Y, Z, zdir='z',  cmap=cm.coolwarm)
    #ax.contour(X, Y, Z, zdir='x',  cmap=cm.coolwarm)
    #ax.contour(X, Y, Z, zdir='y',  cmap=cm.coolwarm)

    #ax.set_xlim3d(1e-6, 1);
    #ax.set_ylim3d(0, 3*pi);
    #ax.set_zlim3d(-pi, 2*pi);
    ax.set_xlabel(r'k$_{h}$ [$\frac{\textrm{Mpc}}{h}$]')
    ax.set_ylabel(r'z')
    ax.set_zlabel(r'T(k,z)');
    
    ax.set_zscale('log'); #ax.xaxis.set_scale('log'); #ax.set_xscale('log') ; #ax.set_yscale('log')
    
    name = out + 'linear_cdm_transfer.png'
    plt.tight_layout()
    plt.savefig(name)
    print(name, 'PRINTED!')

    plt.clf()
    plt.xticks(fontsize=fontsize); plt.yticks(fontsize=fontsize)
    plt.plot(kh, delta[0], linewidth=linewidth)
    plt.xlabel(r'k$_{h}$ [$\frac{Mpc}{h}$]', fontsize=fontsize)
    plt.ylabel(r'T(k,0)', fontsize=24)
    plt.xscale('log') ; plt.yscale('log')
    plt.locator_params(axis='x', numticks=4)
    #plt.xlim([0, 400]);  #plt.ylim([-260, 260])
    name = out + 'linear_cdm_transfer_z0.png'
    plt.tight_layout()
    plt.savefig(name)
    print(name, 'PRINTED!')


def main():
    args = parse_args()
    out = os.path.expanduser(args.out + 'plots/')
    try:
        if not os.path.isdir(out):
            os.makedirs(out)
    except OSError as e:
        print "Ignore OSError from makedirs(work):"
        print e
        pass
    
    cambpath =  os.path.join(args.work)
    camb(cambpath, out)
    
if __name__ == "__main__":
    main()