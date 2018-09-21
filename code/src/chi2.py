def modelvectors(rhos, alpha=0, beta=0, gamma=0, gflag=True, betaflag=True):
    if(gflag and betaflag):
        mvec0 = alpha*rhos[0] + beta*rhos[2] + gamma*rhos[5] 
        mvec1 = alpha*rhos[2] + beta*rhos[1] + gamma*rhos[4] 
        mvec2 = alpha*rhos[5] + beta*rhos[4] + gamma*rhos[3]
        return mvec0, mvec1, mvec2
    elif(gflag and (not betaflag)):
        mvec0 = alpha*rhos[0] + gamma*rhos[5] 
        mvec1 = alpha*rhos[5] + gamma*rhos[3]
        return mvec0, mvec1
    elif((not gflag) and betaflag):
        mvec0 = alpha*rhos[0] + beta*rhos[2]
        mvec1 = alpha*rhos[2] + beta*rhos[1] 
        return mvec0, mvec1
    else:
        mvec0 = alpha*rhos[0]
        return mvec0
        

def modelvar(msigs, gflag=True, betaflag=True):
    #variances sigs^2
    if(gflag and betaflag):
        mvar0 = msigs[0]**2 +msigs[2]**2 + msigs[5]**2 
        mvar1 = msigs[2]**2 +msigs[1]**2 + msigs[4]**2  
        mvar2 = msigs[5]**2 +msigs[4]**2 + msigs[3]**2  
        return mvar0, mvar1, mvar2
    elif(gflag and (not betaflag)):
        mvar0 = msigs[0]**2 +msigs[5]**2
        mvar1 = msigs[5]**2 +msigs[3]**2
        return mvar0,  mvar1
    elif((not gflag) and betaflag):
        mvar0 = msigs[0]**2 +msigs[2]**2
        mvar1 = msigs[2]**2 +msigs[1]**2
        return mvar0,  mvar1
    else:
        mvar0 = msigs[0]**2 
        return mvar0
        
def chi2(model, data,  varmodel, vardata ):
    import numpy as np
    chisq_vec = np.power((model - data), 2)/(varmodel + vardata)
    return chisq_vec.sum()

def CHI2(params, data, eq=None,  gflag=True,  betaflag=True):
    rhos = data['rhos']
    sigrhos = data['sigrhos']
    taus =  data['taus']
    sigtaus = data['sigtaus']
    if(gflag and betaflag):
        #print("Using alpha, beta and gamma")
        alpha, beta, gamma = params
        mvect0, mvect1, mvect2 =  modelvectors(rhos, alpha, beta, gamma, gflag=gflag, betaflag=betaflag)
        mvar0, mvar1, mvar2 = modelvar(sigrhos, gflag=gflag, betaflag=betaflag)
    elif(gflag and (not betaflag)):
        #print("Using alpha and gamma")
        alpha, gamma = params
        mvect0, mvect1=  modelvectors(rhos, alpha=alpha, gamma=gamma, gflag=gflag, betaflag=betaflag)
        mvar0, mvar1= modelvar(sigrhos, gflag=gflag, betaflag=betaflag)
    elif((not gflag) and betaflag):
        #print("Using alpha and beta")
        alpha, beta = params
        mvect0, mvect1=  modelvectors(rhos, alpha=alpha, beta=beta, gflag=gflag, betaflag=betaflag)
        mvar0, mvar1= modelvar(sigrhos, gflag=gflag, betaflag=betaflag)
    else:
        #print("Using only alpha")
        alpha = params
        mvect0 =  modelvectors(rhos, alpha=alpha, gflag=gflag, betaflag=betaflag)
        mvar0 = modelvar(sigrhos, gflag=gflag, betaflag=betaflag)

    dvect0, dvect1, dvect2 = taus[0], taus[1],  taus[2]
    dvar0, dvar1, dvar2 =  sigtaus[0]**2, sigtaus[1]**2, sigtaus[2]**2

    if(not gflag and not betaflag):
        #print("Using all the system equations")
        val=chi2(mvect0, dvect0, mvar0, dvar0 )
        return val
    elif(gflag^betaflag):
        if(eq==0):
            #print("Using first equation")
            val=chi2(mvect0, dvect0, mvar0, dvar0 )
            return val
        elif(eq==1):
            #print("Using second equation")
            val=chi2(mvect1, dvect1, mvar1, dvar1 )
            return val
        else:
            #print("Using all the system of two equations")
            val=chi2(mvect0 + mvect1, dvect0 + dvect1, mvar0 + mvar1, dvar0 + dvar1)
            return val
    else:
        if(eq==0):
            #print("Using first equation")
            val=chi2(mvect0, dvect0, mvar0, dvar0 )
            return val
        elif(eq==1):
            #print("Using second equation")
            val=chi2(mvect1, dvect1, mvar1, dvar1 )
            return val
        elif(eq==2):
            #print("Using third equation")
            val=chi2(mvect2, dvect2, mvar2, dvar2 )
            return val
        else:
            #print("Using all the system of three equations")
            val=chi2(mvect0 + mvect1 + mvect2, dvect0 + dvect1 + dvect2, mvar0 + mvar1 + mvar2, dvar0 + dvar1 + dvar2 )
            return val


def minimize(data, initial_guess, eq=None,  gflag=True, betaflag = True):
    import scipy.optimize as optimize
    result = optimize.minimize(CHI2, initial_guess,args=(data, eq,gflag,betaflag), method='Nelder-Mead', tol=1e-6)
    if result.success:
        fitted_params = result.x
        return fitted_params, result.fun
    else:
        raise ValueError(result.message)
