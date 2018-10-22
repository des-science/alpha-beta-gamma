def modelvector(rhos, params, eq=None, gflag=True, bflag=True):
    if(gflag and bflag):
        alpha, beta, gamma = params
        mvec0 = alpha*rhos[0] + beta*rhos[2] + gamma*rhos[5] 
        mvec1 = alpha*rhos[2] + beta*rhos[1] + gamma*rhos[4] 
        mvec2 = alpha*rhos[5] + beta*rhos[4] + gamma*rhos[3]    
    elif(gflag and (not bflag)):
        alpha, gamma = params
        mvec0 = alpha*rhos[0] + gamma*rhos[5] 
        mvec1 = alpha*rhos[2] + gamma*rhos[4]
        mvec2 = alpha*rhos[5] + gamma*rhos[3]
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
    else:
        return mvec0 +  mvec1 +  mvec2
        
def modelvar(msigs, params, eq =None, gflag=True, bflag=True):
    #variances sigs^2
    if(gflag and bflag):
        alpha, beta, eta = params
        mvar0 = (alpha**2)*msigs[0]**2 +(beta**2)*msigs[2]**2 + (eta**2)*msigs[5]**2 
        mvar1 = (alpha**2)*msigs[2]**2 +(beta**2)*msigs[1]**2 + (eta**2)*msigs[4]**2  
        mvar2 = (alpha**2)*msigs[5]**2 +(beta**2)*msigs[4]**2 + (eta**2)*msigs[3]**2  
    elif(gflag and (not bflag)):
        alpha, eta = params
        mvar0 = (alpha**2)*msigs[0]**2 +(eta**2)*msigs[5]**2
        mvar1 = (alpha**2)*msigs[2]**2 +(eta**2)*msigs[4]**2
        mvar2 = (alpha**2)*msigs[5]**2 +(eta**2)*msigs[3]**2   
    elif((not gflag) and bflag):
        alpha, beta = params
        mvar0 = (alpha**2)*msigs[0]**2 +(beta**2)*msigs[2]**2
        mvar1 = (alpha**2)*msigs[2]**2 +(beta**2)*msigs[1]**2
        mvar2 = (alpha**2)*msigs[5]**2 +(beta**2)*msigs[4]**2
    else:
        alpha = params
        mvar0 = (alpha**2)*msigs[0]**2
        mvar1 = (alpha**2)*msigs[2]**2
        mvar2 = (alpha**2)*msigs[5]**2 
    if(eq==0):
        return mvar0
    elif(eq==1):
        return mvar1
    elif(eq==2):
        return mvar2
    else:
        return mvar0 +  mvar1 +  mvar2
   
def datavector(taus, eq=None):
    if(eq==0):
        return taus[0]
    elif(eq==1):
        return taus[1]
    elif(eq==2):
        return taus[2]
    else:
        return taus[0] + taus[1] + taus[2]
def datavar(sigtaus, eq=None):
    if(eq==0):
        return sigtaus[0]**2
    elif(eq==1):
        return sigtaus[1]**2
    elif(eq==2):
        return sigtaus[2]**2
    else:
        return sigtaus[0]**2 + sigtaus[1]**2 + sigtaus[2]**2
    
def chi2(modelvec, datavec,  varmodel, vardata,  moderr=False ):
    import numpy as np
    if(moderr):
        chisq_vec = np.power((modelvec - datavec), 2)/(varmodel + vardata)
    else:
        chisq_vec = np.power((modelvec - datavec), 2)/(vardata)
    return chisq_vec.sum()

def CHI2(params, data, eq=None,  gflag=True,  bflag=True, moderr=False):
    rhos = data['rhos'];sigrhos = data['sigrhos']
    taus =  data['taus'];sigtaus = data['sigtaus']
    dvect=  datavector(taus, eq=eq)
    dvar = datavar(sigtaus, eq=eq)
    mvect=  modelvector(rhos, params, eq=eq,  gflag=gflag, bflag=bflag)
    mvar = modelvar(sigrhos, params, eq=eq, gflag=gflag, bflag=bflag)
    val=chi2(mvect, dvect, mvar, dvar, moderr=moderr )
    return val
  
def CHI2shifted(params,data, svalue, eq=None,  gflag=True,  bflag=True,  moderr=False):
    return CHI2(params, data, eq=eq, gflag=gflag,bflag=bflag, moderr=moderr) -svalue 

def plotCHI2(pars, data,x_arr, filename, eq=None, gflag=True,  bflag=True ):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    import numpy as np
    if(gflag and bflag):
        print("")
    elif(gflag and (not bflag)):
        print("")
    elif((not gflag) and bflag):
        print("")
    else:
        xmin, xmax, npoints = x_arr
        alphas = np.linspace(xmin, xmax, npoints)
        chisq = [CHI2(x,data,eq=eq, gflag=False, bflag=False ) for x in alphas]
        plt.clf()
        plt.xlabel(r"$\alpha$")
        plt.ylabel(r"$\chi^{2}$")
        plt.xlim([xmin, xmax])
        plt.plot(alphas, chisq)
        print("Printing file: ", filename)
        plt.savefig(filename)#, dpi=150)
        
def plotCHI2shifted(pars, data, x_arr, svalue, filename , eq=None, gflag=True,  bflag=True):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    import numpy as np
    if(gflag and bflag):
        print("")
    elif(gflag and (not bflag)):
        print("")
    elif((not gflag) and bflag):
        print("")
    else:
        xmin, xmax, npoints = x_arr
        alphas = np.linspace(xmin, xmax, npoints)
        chisq = [CHI2shifted(x,data,svalue, eq=eq, gflag=False, bflag=False ) for x in alphas]
        plt.clf()
        plt.xlabel(r"$\alpha$")
        plt.ylabel(r"$\chi^{2}$")
        plt.xlim([xmin, xmax])
        plt.plot(alphas, chisq)
        print("Printing file: ",  filename)
        plt.savefig(filename)#, dpi=150)
        
def minimizeCHI2(data, initial_guess, eq=None,  gflag=True, bflag = True, moderr=False):
    import scipy.optimize as optimize
    result = optimize.minimize(CHI2, initial_guess,args=(data,eq,gflag,bflag, moderr), method='Nelder-Mead', tol=1e-6)
    if result.success:
        fitted_params = result.x
        return fitted_params, result.fun
    else:
        raise ValueError(result.message)
