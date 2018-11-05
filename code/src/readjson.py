import json

def read_rhos(stat_file, maxscale=None,  maxbin =None ):
    import numpy as np

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

def read_taus(stat_file, maxscale=None, maxbin =None ):
    import numpy as np

    # Read the json file 
    with open(stat_file,'r') as f:
        stats = json.load(f)
        
    #print' stats = ',stats
    if len(stats) == 1:  # I used to save a list of length 1 that in turn was a list
        stats = stats[0]
            
        ( meanlogr,
          tau0p,
          tau0p_im,
          tau0m,
          tau0m_im,
          var0,
          tau2p,
          tau2p_im,
          tau2m,
          tau2m_im,
          var2,
          tau5p,
          tau5p_im,
          tau5m,
          tau5m_im,
          var5,
        ) = stats[:16]

        #Finally this are the arrays with the data
        meanr = np.exp(meanlogr)
        tau0p = np.array(tau0p)
        tau0m = np.array(tau0m)
        tau2p = np.array(tau2p)
        tau2m = np.array(tau2m)
        tau5p = np.array(tau5p)
        tau5m = np.array(tau5m)
        sig_tau0 = np.sqrt(var0)
        sig_tau2 = np.sqrt(var2)
        sig_tau5 = np.sqrt(var5)

        if(maxscale):
            maxs = maxscale
            meanr = meanr[meanr<maxs]
            idx = len(meanr)
            tau0p = tau0p[:idx]; tau2p = tau2p[:idx]; tau5p = tau5p[:idx];
            sig_tau0 = sig_tau0[:idx]; sig_tau2 = sig_tau2[:idx]; sig_tau5 = sig_tau5[:idx]

        if(maxbin):
            meanr = meanr[: maxbin]; tau0p = tau0p[:maxbin];
            tau2p = tau2p[:maxbin]; tau5p = tau5p[:maxbin];
            sig_tau0 = sig_tau0[:maxbin]; sig_tau2 = sig_tau2[:maxbin]; sig_tau5 = sig_tau5[:maxbin]

        return meanr, tau0p, tau2p, tau5p, sig_tau0, sig_tau2, sig_tau5 
              

 
 
