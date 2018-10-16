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
    parser.add_argument('--bands', default='grizY', type=str,
                         help='Limit to the given bands')
  
   
    
    
    args = parser.parse_args()

    return args

def main():
    import numpy as np
    from read_psf_cats import read_data,  toList
    from readjson import read_rhos

    args = parse_args()
    #usualrhos =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/used/rho_all_reserved_irz_usual.json"
    usualrhos = "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_unmod_irz.json"
    #modrhos = "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/used/rho_all_reserved_irz_modified.json"
    modrhos =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/rho_all_reserved_mod_irz.json"
    umeanr, urho0p, urho1p, urho2p, urho3p, urho4p, urho5p, usig_rho0, usig_rho1, usig_rho2, usig_rho3, usig_rho4, usig_rho5 = read_rhos(usualrhos)
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 = read_rhos(modrhos)

    print("<ep'ep'> - <epep>")
    print(rho0p - urho0p)
    print("<q'q'> - <qq>")
    print(rho1p - urho1p)
    print("<w'w'>-<ww>")
    print(rho3p - urho3p)

    keys = ['ra', 'dec','obs_e1', 'obs_e2', 'obs_T',
            'piff_e1', 'piff_e2', 'piff_T']
 
    exps = toList(args.exps_file)
    data, bands, tilings = read_data(exps, args.piff_cat , keys,
                                     limit_bands=args.bands,
                                     use_reserved=True)

    prefix = 'piff'
    e1 = data['obs_e1']
    e2 = data['obs_e2']
    
    p_e1 = data[prefix+'_e1']
    p_e2 = data[prefix+'_e2']
    
    de1 = e1-p_e1
    de2 = e2-p_e2

    
    T = data['obs_T']
    p_T = data[prefix+'_T']
    dt = (T-p_T)/T
    w1 = p_e1*dt
    w2 = p_e2*dt


    print("<p_e>^2")
    meanpe = np.mean(p_e1)**2 + np.mean(p_e2)**2
    meanpe2 = np.mean(np.sqrt(p_e1**2 + p_e2**2))**2
    print(meanpe, meanpe2)
    print("<de>^2")
    meande = np.mean(de1)**2 + np.mean(de2)**2
    meande2 = np.mean(np.sqrt(de1**2 + de2**2))**2
    print(meande, meande2 )
    print("<w>^2")
    meanw = np.mean(w1)**2 + np.mean(w2)**2
    meanw2 = np.mean(np.sqrt(w1**2 + w2**2))**2
    print(meanw, meanw2)
    
    
    
if __name__ == "__main__":
    main()
