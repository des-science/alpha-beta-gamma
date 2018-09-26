def modelvectors(rhos, alpha=0, beta=0, gamma=0, gflag=True, bflag=True):
    if(gflag and bflag):
        mvec0 = alpha*rhos[0] + beta*rhos[2] + gamma*rhos[5] 
        mvec1 = alpha*rhos[2] + beta*rhos[1] + gamma*rhos[4] 
        mvec2 = alpha*rhos[5] + beta*rhos[4] + gamma*rhos[3]
        return mvec0, mvec1, mvec2
    elif(gflag and (not bflag)):
        mvec0 = alpha*rhos[0] + gamma*rhos[5] 
        mvec1 = alpha*rhos[5] + gamma*rhos[3]
        return mvec0, mvec1
    elif((not gflag) and bflag):
        mvec0 = alpha*rhos[0] + beta*rhos[2]
        mvec1 = alpha*rhos[2] + beta*rhos[1] 
        return mvec0, mvec1
    else:
        mvec0 = alpha*rhos[0]
        return mvec0
        

def modelvar(msigs, gflag=True, bflag=True):
    #variances sigs^2
    if(gflag and bflag):
        mvar0 = msigs[0]**2 +msigs[2]**2 + msigs[5]**2 
        mvar1 = msigs[2]**2 +msigs[1]**2 + msigs[4]**2  
        mvar2 = msigs[5]**2 +msigs[4]**2 + msigs[3]**2  
        return mvar0, mvar1, mvar2
    elif(gflag and (not bflag)):
        mvar0 = msigs[0]**2 +msigs[5]**2
        mvar1 = msigs[5]**2 +msigs[3]**2
        return mvar0,  mvar1
    elif((not gflag) and bflag):
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

def CHI2(params, data, eq=None,  gflag=True,  bflag=True):
    rhos = data['rhos']
    sigrhos = data['sigrhos']
    taus =  data['taus']
    sigtaus = data['sigtaus']
    if(gflag and bflag):
        #print("Using alpha, beta and gamma")
        alpha, beta, gamma = params
        mvect0, mvect1, mvect2 =  modelvectors(rhos, alpha, beta, gamma, gflag=gflag, bflag=bflag)
        mvar0, mvar1, mvar2 = modelvar(sigrhos, gflag=gflag, bflag=bflag)
    elif(gflag and (not bflag)):
        #print("Using alpha and gamma")
        alpha, gamma = params
        mvect0, mvect1=  modelvectors(rhos, alpha=alpha, gamma=gamma, gflag=gflag, bflag=bflag)
        mvar0, mvar1= modelvar(sigrhos, gflag=gflag, bflag=bflag)
    elif((not gflag) and bflag):
        #print("Using alpha and beta")
        alpha, beta = params
        mvect0, mvect1=  modelvectors(rhos, alpha=alpha, beta=beta, gflag=gflag, bflag=bflag)
        mvar0, mvar1= modelvar(sigrhos, gflag=gflag, bflag=bflag)
    else:
        #print("Using only alpha")
        alpha = params
        mvect0 =  modelvectors(rhos, alpha=alpha, gflag=gflag, bflag=bflag)
        mvar0 = modelvar(sigrhos, gflag=gflag, bflag=bflag)

    dvect0, dvect1, dvect2 = taus[0], taus[1],  taus[2]
    dvar0, dvar1, dvar2 =  sigtaus[0]**2, sigtaus[1]**2, sigtaus[2]**2

    if(not gflag and not bflag):
        #print("Using all the system equations")
        val=chi2(mvect0, dvect0, mvar0, dvar0 )
        return val
    elif(gflag^bflag):
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

def CHI2shifted(params,data, svalue, eq=None,  gflag=True,  bflag=True):
    return CHI2(params, data, eq=eq, gflag=gflag,bflag=bflag) -svalue 

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
        
def minimizeCHI2(data, initial_guess, eq=None,  gflag=True, bflag = True):
    import scipy.optimize as optimize
    result = optimize.minimize(CHI2, initial_guess,args=(data,eq,gflag,bflag), method='Nelder-Mead', tol=1e-6)
    if result.success:
        fitted_params = result.x
        return fitted_params, result.fun
    else:
        raise ValueError(result.message)
