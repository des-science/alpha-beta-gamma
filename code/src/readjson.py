import json

def read_rhos(stat_file ):
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
        return meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 

def read_sigmas(stat_file ):
    import numpy as np

    # Read the json file 
    with open(stat_file,'r') as f:
        stats = json.load(f)
        
    #print' stats = ',stats
    if len(stats) == 1:  # I used to save a list of length 1 that in turn was a list
        stats = stats[0]
            
        ( meanlogr,
          sigma0p,
          sigma0p_im,
          sigma0m,
          sigma0m_im,
          var0,
          sigma2p,
          sigma2p_im,
          sigma2m,
          sigma2m_im,
          var2,
          sigma5p,
          sigma5p_im,
          sigma5m,
          sigma5m_im,
          var5,
        ) = stats[:16]

        #Finally this are the arrays with the data
        meanr = np.exp(meanlogr)
        sigma0p = np.array(sigma0p)
        sigma0m = np.array(sigma0m)
        sigma2p = np.array(sigma2p)
        sigma2m = np.array(sigma2m)
        sigma5p = np.array(sigma5p)
        sigma5m = np.array(sigma5m)
        sig_sigma0 = np.sqrt(var0)
        sig_sigma2 = np.sqrt(var2)
        sig_sigma5 = np.sqrt(var5)
        return meanr, sigma0p, sigma2p, sigma5p, sig_sigma0, sig_sigma2, sig_sigma5 
              

 
 
