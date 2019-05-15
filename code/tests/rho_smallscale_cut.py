#run abn test and save files arrays of parameters going out to a maximum bin. 
#plotting each term in the equation of correlations, to see if there are reason to cancel some part of the model.
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('SVA1StyleSheet.mplstyle')

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Correlation of reserved stars')
    
    parser.add_argument('--piff_cat',
                        default='/home2/dfa/sobreira/alsina/catalogs/y3a1-v29',
                        help='Full Path to the Only stars Piff catalog')
    parser.add_argument('--exps_file',
                        default='/home/dfa/sobreira/alsina/DESWL/psf/ally3.grizY',
                        #default='/home/dfa/sobreira/alsina/DESWL/psf/testexp', 
                        help='list of exposures (in lieu of separate exps)')
    parser.add_argument('--bands', default='riz', type=str,
                         help='Limit to the given bands')
    parser.add_argument('--bandcombo', default=False,
                        action='store_const', const=True,
                        help='run rho2 for all combination of bands, if false run particular combination defined in band')
    parser.add_argument('--use_reserved', default=True,
                        action='store_const', const=True,
                        help='just use the objects with the RESERVED flag')
    parser.add_argument('--frac', default=1., type=float,
                        help='Choose a random fraction of the input stars')
    parser.add_argument('--mod', default=True,
                        action='store_const', const=True,
                        help='If true it substracts the mean to each field before calculate correlations')
    parser.add_argument('--sn', default=True,
                        action='store_const', const=True,
                        help='If true multiply by 2 the variances of all correlations. Shape-noise error.')
    parser.add_argument('--obs', default=False,
                        action='store_const', const=True,
                        help='Use e_obs instead of e_piff to calculate modified rho stats')
    parser.add_argument('--outpath', default='/home/dfa/sobreira/alsina/alpha-beta-gamma/code/tests/rhominscale/',
                        help='location of the output of the files')
    args = parser.parse_args()
    return args

def write_stats( stat_file, rho0, rho1, rho2, rho3, rho4, rho5):
    import json
    stats = [
        rho0.meanlogr.tolist(),
        rho0.xip.tolist(),
        rho0.xip_im.tolist(),
        rho0.xim.tolist(),
        rho0.xim_im.tolist(),
        rho0.varxi.tolist(),
        rho1.xip.tolist(),
        rho1.xip_im.tolist(),
        rho1.xim.tolist(),
        rho1.xim_im.tolist(),
        rho1.varxi.tolist(),
        rho2.xip.tolist(),
        rho2.xip_im.tolist(),
        rho2.xim.tolist(),
        rho2.xim_im.tolist(),
        rho2.varxi.tolist(),
        rho3.xip.tolist(),
        rho3.xip_im.tolist(),
        rho3.xim.tolist(),
        rho3.xim_im.tolist(),
        rho3.varxi.tolist(),
        rho4.xip.tolist(),
        rho4.xip_im.tolist(),
        rho4.xim.tolist(),
        rho4.xim_im.tolist(),
        rho4.varxi.tolist(),
        rho5.xip.tolist(),
        rho5.xip_im.tolist(),
        rho5.xim.tolist(),
        rho5.xim_im.tolist(),
        rho5.varxi.tolist(),
    ]
    #print('stats = ',stats)
    print('stat_file = ',stat_file)
    with open(stat_file,'w') as fp:
        json.dump([stats], fp)
    print('Done writing ',stat_file)

def read_rhos(stat_file, maxscale=None,  maxbin =None ):
    # Read the json file 
    with open(stat_file,'r') as f:
        stats = json.load(f)
        
    #print' stats = ',stats
    if len(stats) == 1:  # I used to save a list of length 1 that in turn was a list
        stats = stats[0]
            
        ( meanlogr,
          rho0p,
          rho0p_im,
          rho0m,
          rho0m_im,
          var0,
          rho1p,
          rho1p_im,
          rho1m,
          rho1m_im,
          var1,
          rho2p,
          rho2p_im,
          rho2m,
          rho2m_im,
          var2,
          rho3p,
          rho3p_im,
          rho3m,
          rho3m_im,
          var3,
          rho4p,
          rho4p_im,
          rho4m,
          rho4m_im,
          var4,
          rho5p,
          rho5p_im,
          rho5m,
          rho5m_im,
          var5,
        ) = stats[:31]

        #Finally this are the arrays with the data
        meanr = np.exp(meanlogr)
        rho0p = np.array(rho0p)
        rho0m = np.array(rho0m)
        rho1p = np.array(rho1p)
        rho1m = np.array(rho1m)
        rho2p = np.array(rho2p)
        rho2m = np.array(rho2m)
        rho3p = np.array(rho3p)
        rho3m = np.array(rho3m)
        rho4p = np.array(rho4p)
        rho4m = np.array(rho4m)
        rho5p = np.array(rho5p)
        rho5m = np.array(rho5m)
        sig_rho0 = np.sqrt(var0)
        sig_rho1 = np.sqrt(var1)
        sig_rho2 = np.sqrt(var2)
        sig_rho3 = np.sqrt(var3)
        sig_rho4 = np.sqrt(var4)
        sig_rho5 = np.sqrt(var5)

        if(maxscale):
            maxs = maxscale
            meanr = meanr[meanr<maxs]
            idx = len(meanr)
            rho0p = rho0p[:idx]; rho1p = rho1p[:idx]; rho2p = rho2p[:idx]; rho3p = rho3p[:idx];
            rho4p = rho4p[:idx]; rho5p = rho5p[:idx]; sig_rho0 = sig_rho0[:idx];
            sig_rho1 = sig_rho1[:idx]; sig_rho2 = sig_rho2[:idx]; sig_rho3 = sig_rho3[:idx];
            sig_rho4 = sig_rho4[:idx]; sig_rho5 = sig_rho5[:idx]

        if(maxbin):
            meanr = meanr[:maxbin]
            rho0p = rho0p[:maxbin]; rho1p = rho1p[:maxbin]; rho2p = rho2p[:maxbin];
            rho3p = rho3p[:maxbin]; rho4p = rho4p[:maxbin]; rho5p = rho5p[:maxbin];
            sig_rho0 = sig_rho0[:maxbin]; sig_rho1 = sig_rho1[:maxbin];
            sig_rho2 = sig_rho2[:maxbin]; sig_rho3 = sig_rho3[:maxbin];
            sig_rho4 = sig_rho4[:maxbin]; sig_rho5 = sig_rho5[:maxbin]
        
        return meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5     
def measure_rho(data, min_sep = 0.5,  max_sep=300, sep_units='arcmin',  tag=None, prefix='piff', mod=True,  obs=False ):
    """Compute the rho statistics
    """
    import treecorr
    import numpy as np

    e1 = data['obs_e1']
    e2 = data['obs_e2']
    p_e1 = data[prefix+'_e1']
    p_e2 = data[prefix+'_e2']
    T = data['obs_T']
    p_T = data[prefix+'_T']

    de1 = e1-p_e1
    de2 = e2-p_e2
    dt = (T-p_T)/T
    w1 = p_e1*dt
    w2 = p_e2*dt
    w1obs = e1*dt 
    w2obs = e2*dt 
        
    
    #Modified ellipticities
    if(mod):
        e1 = e1 - np.array(np.mean(e1))
        e2 = e2 - np.array(np.mean(e2))       
        p_e1 = p_e1 - np.array(np.mean(p_e1))
        p_e2 = p_e2 - np.array(np.mean(p_e2))
        de1 = de1 - np.array(np.mean(de1))
        de2 = de2 - np.array(np.mean(de2))
        w1 = w1 - np.array(np.mean(w1))
        w2 = w2 - np.array(np.mean(w2))
        w1obs = w1obs - np.array(np.mean(w1obs))
        w2obs = w2obs - np.array(np.mean(w2obs))
        
    ra = data['ra']
    dec = data['dec']
    print('ra = ',ra)
    print('dec = ',dec)
    if(obs):
        ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=e1, g2=e2)
        decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
        wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1obs, g2=w2obs)
    else:
        ecat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=p_e1, g2=p_e2)
        decat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=de1, g2=de2)
        #wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1, g2=w2)
        wcat = treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', g1=w1obs, g2=w2obs)
    ecat.name = 'ecat'
    decat.name = 'decat'
    wcat.name = 'wcat'
    if tag is not None:
        for cat in [ ecat, decat, wcat ]:
            cat.name = tag + ":"  + cat.name

    #bin_config = dict(sep_units = sep_units, nbins = 20, min_sep = 2.5, max_sep = 250,)
    bin_config = dict(sep_units = sep_units, bin_slop = 0.1, min_sep = min_sep, max_sep = max_sep, bin_size = 0.2)
 
    

    results = []
    for (cat1, cat2) in [(ecat, ecat), 
                         (decat, decat),
                          (decat, ecat),
                          (wcat, wcat),
                          (decat, wcat),
                          (ecat, wcat) ]:
        print('Doing correlation of %s vs %s'%(cat1.name, cat2.name))

        rho = treecorr.GGCorrelation(bin_config, verbose=2)

        if cat1 is cat2:
            rho.process(cat1)
        else:
            rho.process(cat1, cat2)
        print('mean xi+ = ',rho.xip.mean())
        print('mean xi- = ',rho.xim.mean())
        results.append(rho)

    return results

def run_rhos(args,  outpath ):
    from read_psf_cats import read_data,  toList
    #STATISTIC USING ONLY RESERVED STARS
    keys = ['ra', 'dec','obs_e1', 'obs_e2', 'obs_T',
            'piff_e1', 'piff_e2', 'piff_T',  'mag']
 
    exps = toList(args.exps_file)
    data, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=args.use_reserved)
    print("Objects",  len(data))
    data = data[data['mag']<20]
    print("Objects with magnitude <20",  len(data))

    min_sep_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for min_sep in min_sep_list:
        stats = measure_rho(data, min_sep=min_sep, mod=args.mod,  obs=args.obs)
        stat_file = os.path.join(outpath, "rho_%f.json"%(min_sep))
        write_stats(stat_file,*stats)

def plot_correlations(outpath):

    from plot_stats import pretty_rho
    import numpy as np
    
    ylims =  [[ - 0.1, 0.1],[ - 30, 30],[ -600, 600] ] 
    outputnames = ['rho0.png', 'rho1.png', 'rho2.png', 'rho3.png', 'rho4.png', 'rho5.png']
    colors = ['black', 'green', 'blue', 'red', 'gray', 'pink']
    min_sep_list = [0.05, 0.1, 0.2, 0.3, 0.4]
    #min_sep_list = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for min_sep in enumeratemin_sep_list:
        meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 =read_rhos(os.path.join(outpath, "rho_%f.json"%(min_sep)))

        plt.clf()
        plt.figure(0)
        pretty_rho(meanr, rho0p, sig_rho0, legend=str(min_sep), lfontsize=24,
                   color='black', marker='o',
                   ylabel=r'$\rho_{0}(\theta)$',title=r'$\rho_{0}(\theta)$',
                   xlim=None, ylim=None)

        plt.clf()
        plt.figure(1)
        pretty_rho(meanr, rho1p, sig_rho1, legend=str(min_sep), lfontsize=24,
                   color='black', marker='o',
                   ylabel=r'$\rho_{1}(\theta)$',title=r'$\rho_{1}(\theta)$',
                   xlim=None, ylim=None)

        plt.clf()
        plt.figure(2)
        pretty_rho(meanr, rho2p, sig_rho2, legend=str(min_sep), lfontsize=24,
                   color='black', marker='o',
                   ylabel=r'$\rho_{2}(\theta)$',title=r'$\rho_{2}(\theta)$',
                   xlim=None, ylim=None)

        plt.clf()
        plt.figure(3)
        pretty_rho(meanr, rho3p, sig_rho3, legend=str(min_sep), lfontsize=24,
                   color='black', marker='o',
                   ylabel=r'$\rho_{3}(\theta)$',title=r'$\rho_{3}(\theta)$',
                   xlim=None, ylim=None)

        plt.clf()
        plt.figure(4)
        pretty_rho(meanr, rho4p, sig_rho4, legend=str(min_sep), lfontsize=24,
                   color='black', marker='o',
                   ylabel=r'$\rho_{4}(\theta)$',title=r'$\rho_{4}(\theta)$',
                   xlim=None, ylim=None)

        plt.clf()
        plt.figure(5)
        pretty_rho(meanr, rho5p, sig_rho5, legend=str(min_sep), lfontsize=24,
                   color='black', marker='o',
                   ylabel=r'$\rho_{5}(\theta)$',title=r'$\rho_{5}(\theta)$',
                   xlim=None, ylim=None)
        

    for i in range(6):
        plt.figure(i)
        plt.tight_layout()
        print("Printing :", outpath + outputnames[t])
        plt.savefig(outpath +outputnames[i],  dpi=200)
        

       
        
        
def main():
    import sys
    sys.path.insert(0, '/home/dfa/sobreira/alsina/alpha-beta-gamma/code/src')
    #sys.path.insert(0, '/global/cscratch1/sd/alsina/alpha-beta-gamma/code/src')

    

    
    args = parse_args()

    #Make directory where the ouput data will be
    outpath = os.path.expanduser(args.outpath)
    try:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
    except OSError:
        if not os.path.exists(outpath): raise

    run_rhos(args, outpath)
        

    
    


    
if __name__ == "__main__":
    main()
