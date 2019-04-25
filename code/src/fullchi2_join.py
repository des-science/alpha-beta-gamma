def modelvector(rhos, params, eq=None, gflag=True, bflag=True):
    if(gflag and bflag):
        alpha, beta, eta = params
        mvec0 = alpha*rhos[0] + beta*rhos[2] + eta*rhos[5] 
        mvec1 = alpha*rhos[2] + beta*rhos[1] + eta*rhos[4] 
        mvec2 = alpha*rhos[5] + beta*rhos[4] + eta*rhos[3]    
    elif(gflag and (not bflag)):
        alpha, eta = params
        mvec0 = alpha*rhos[0] + eta*rhos[5] 
        mvec1 = alpha*rhos[2] + eta*rhos[4]
        mvec2 = alpha*rhos[5] + eta*rhos[3]
    elif((not gflag) and bflag):
        alpha, beta = params
        mvec0 = alpha*rhos[0] + beta*rhos[2]
        mvec1 = alpha*rhos[2] + beta*rhos[1]
        mvec2 = alpha*rhos[5] + beta*rhos[4] 
    else:
        alpha = params
        mvec0 = alpha*rhos[0]
        mvec1 = alpha*rhos[2]
        mvec2 = alpha*rhos[5]
    if(eq==0):
        return mvec0
    elif(eq==1):
        return mvec1
    elif(eq==2):
        return mvec2
    elif(eq==[0,1]):
        return mvec0 + mvec1
    elif(eq==[0,2]):
        return mvec0 + mvec2
    elif(eq==[1,2]):
        return mvec1 + mvec2    
    else:
        return mvec0 +  mvec1 +  mvec2
        
def modelcov(covrhos, params, eq =None, gflag=True, bflag=True):
    #variances sigs^2
    if(gflag and bflag):
        alpha, beta, eta = params
        mvar0 = (alpha**2)*covrhos[0]+(beta**2)*covrhos[2]+ (eta**2)*covrhos[5]
        mvar1 = (alpha**2)*covrhos[2] +(beta**2)*covrhos[1] + (eta**2)*covrhos[4] 
        mvar2 = (alpha**2)*covrhos[5] +(beta**2)*covrhos[4] + (eta**2)*covrhos[3] 
    elif(gflag and (not bflag)):
        alpha, eta = params
        mvar0 = (alpha**2)*covrhos[0] +(eta**2)*covrhos[5]
        mvar1 = (alpha**2)*covrhos[2] +(eta**2)*covrhos[4]
        mvar2 = (alpha**2)*covrhos[5] +(eta**2)*covrhos[3]  
    elif((not gflag) and bflag):
        alpha, beta = params
        mvar0 = (alpha**2)*covrhos[0] +(beta**2)*covrhos[2]
        mvar1 = (alpha**2)*covrhos[2] +(beta**2)*covrhos[1]
        mvar2 = (alpha**2)*covrhos[5] +(beta**2)*covrhos[4]
    else:
        alpha = params
        mvar0 = (alpha**2)*covrhos[0]
        mvar1 = (alpha**2)*covrhos[2]
        mvar2 = (alpha**2)*covrhos[5]
    if(eq==0):
        return mvar0
    elif(eq==1):
        return mvar1
    elif(eq==2):
        return mvar2
    elif(eq==[0,1]):
        return mvar0 + mvar1
    elif(eq==[0,2]):
        return mvar0 + mvar2
    elif(eq==[1,2]):
        return mvar1 + mvar2  
    else:
        return mvar0 +  mvar1 +  mvar2
   
def datavector(taus, eq=None):
    if(eq==0):
        return taus[0]
    elif(eq==1):
        return taus[1]
    elif(eq==2):
        return taus[2]
    elif(eq==[0,1]):
        return taus[0] + taus[1]
    elif(eq==[0,2]):
        return taus[0] + taus[2]
    elif(eq==[1,2]):
        return taus[1] + taus[2]
    else:
        return taus[0] + taus[1] + taus[2]
def datacov(covtaus, eq=None):
    if(eq==0):
        return covtaus[0]
    elif(eq==1):
        return covtaus[1]
    elif(eq==2):
        return covtaus[2]
    elif(eq==[0,1]):
        return covtaus[0] + covtaus[1]
    elif(eq==[0,2]):
        return covtaus[0] + covtaus[2]
    elif(eq==[1,2]):
        return covtaus[1] + covtaus[2]
    else:
        return covtaus[0] + covtaus[1] + covtaus[2]
     
def CHI2(params, data, eq=None,  gflag=True,  bflag=True, moderr=False):
    rhosp = data['rhosp'];covrhosp = data['covrhosp']
    rhosm = data['rhosm'];covrhosm = data['covrhosm']
    tausp =  data['tausp'];covtausp = data['covtausp']
    tausm =  data['tausm'];covtausm = data['covtausm']
    dvect=  datavector(tausp, eq=eq)
    dvect+= datavector(tausm, eq=eq)
    mvect=  modelvector(rhosp,params,eq=eq,gflag=gflag,bflag=bflag)
    mvect+= modelvector(rhosm,params,eq=eq,gflag=gflag,bflag=bflag)
    dcov_mat = datacov(covtausp, eq=eq)
    dcov_mat+= datacov(covtausm, eq=eq)
    mcov_mat = modelcov(covrhosp, params, eq=eq, gflag=gflag, bflag=bflag)
    mcov_mat+= modelcov(covrhosm, params, eq=eq, gflag=gflag, bflag=bflag)
    val=chi2(mvect, dvect, mcov_mat, dcov_mat, moderr=moderr )
    return val
def chi2(modelvec, datavec,  covmodel, covdata,  moderr=False ):
    import numpy as np
    d =  np.array([modelvec - datavec])
    if(moderr):
        cov_inv = np.linalg.inv(covdata + covmodel)
    else:
        cov_inv = np.linalg.inv(covdata)
        
    chisq = np.dot(np.dot(d,cov_inv), d.T)
    return chisq[0][0]
    
def minimizeCHI2(data, initial_guess, eq=None,  gflag=True, bflag = True, moderr=False):
    import scipy.optimize as optimize
    result = optimize.minimize(CHI2, initial_guess,args=(data,eq,gflag,bflag, moderr), method='Nelder-Mead', tol=1e-6)
    if result.success:
        fitted_params = result.x
        return fitted_params, result.fun
    else:
        raise ValueError(result.message)
