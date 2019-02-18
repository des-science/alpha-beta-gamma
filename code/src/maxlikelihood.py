#Log natural of the prior 
def logprior(pars, gflag=True,bflag = True):
    import numpy as np
    l = 1000
    alpha_min, beta_min, eta_min  = -l,-l,-l
    alpha_max, beta_max, eta_max  = l,l,l
    if(gflag and bflag):
        #print("Using alpha, beta and delta")
        alpha, beta, eta = pars
        if ( alpha_min<alpha< alpha_max)and( beta_min < beta< beta_max)and( eta_min<eta<eta_max):
            return 0.0
        return -np.inf
    elif(gflag and (not bflag)):
        #print("Using alpha and gamma")
        alpha, eta = pars
    elif((not gflag) and bflag):
        #print("Using alpha and beta")
        alpha, beta = pars
        if (alpha_min<alpha<alpha_max)and( beta_min< beta< beta_max):
            return 0.0
        return -np.inf
    else:
        #print("Using only alpha")
        alpha = pars
        if ( alpha_min < alpha < alpha_max):
            return 0.0
        return -np.inf
##Log natural of the likelihood function. Gaussian.
def loglike(chisq):
    return -0.5*chisq
##Log natural of the posterior
def logpost(pars, data,svalue=None, eq=None, gflag=True,bflag = True, moderr=False):
    import numpy as np
    from chi2 import CHI2, CHI2shifted
    if(svalue):
        chisq = CHI2shifted(pars, data,svalue, eq=eq, gflag=gflag, bflag=bflag, moderr=moderr )
    else:
        chisq = CHI2(pars, data,eq=eq, gflag=gflag, bflag=bflag, moderr=moderr )
    lp = logprior(pars, gflag=gflag,bflag = bflag)
    if not np.isfinite(lp):
        return -np.inf
    return lp + loglike(chisq)

def alpha_percentil(p, mu, data, svalue=None, eq=None):
    from scipy.integrate import quad
    from scipy.optimize import brentq
    import numpy as np
    max_loglike= logpost(mu,data, svalue=svalue, eq=None, gflag=False, bflag=False, moderr=False)
    print(mu,  max_loglike)
    def cumulative_like(alpha):
        integral, _ = quad(lambda x: np.exp(logpost(x, data, svalue=svalue,  eq=eq, gflag=False, bflag=False, moderr=False ) - max_loglike), -np.inf, alpha)
        return integral
    total_int = cumulative_like(np.inf)
    #print(total_int)
    perc = brentq(lambda x: cumulative_like(x) - p*total_int, -np.inf, np.inf)
    return perc

def corner_plot(samples, labels, title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    import corner
    import numpy as np
    burn = 5000
    samples_burned = np.c_[[par[burn:] for par in samples]]
    fig = corner.corner(samples_burned.T, labels=labels,
                        quantiles=[0.16, 0.5, 0.84],  #-1sigma,0sigma,1sigma
                        levels=(1-np.exp(-0.5), 1-np.exp(-2), 1-np.exp(-9./2)), #1sigma, 2sigma and 3sigma contours
                        show_titles=True, title_kwargs={"fontsize": 12},
                        smooth1d=None, plot_contours=True,
                        no_fill_contours=False, plot_density=True,)
    print("Printing file:",  title)
    plt.savefig(title)
    print(title, "Printed")
def MCMC(best_pars,data, nwalkers=50, nsteps=1000, namemc='mcmc.png', namecont='contcurve.png', svalue=None,  eq=None, gflag=True,bflag = True , moderr=False,  plot=True, nsig=1):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    import emcee  
    import numpy as np
    if(gflag and bflag):
        # initial position at maximum likelihood values
        ndim = 3
        pos = [best_pars + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        # MCMC chain with 50 walkers and 1000 steps
        sampler = emcee.EnsembleSampler(nwalkers, ndim, logpost, threads=4,  args=(data,svalue,eq,gflag,bflag, moderr) )
        print("Runing MCMC ...")
        sampler.run_mcmc(pos, nsteps)
        print("Run finished")
        # Getting chains
        alpha_chain = sampler.chain[:,:,0]
        beta_chain = sampler.chain[:,:,1]
        delta_chain = sampler.chain[:,:,2]

        # Reshaping
        alpha_chain_flat = np.reshape(alpha_chain, (nwalkers*nsteps,))
        beta_chain_flat = np.reshape(beta_chain, (nwalkers*nsteps,))
        delta_chain_flat = np.reshape(delta_chain, (nwalkers*nsteps,))

        samples = np.c_[alpha_chain_flat, beta_chain_flat, delta_chain_flat].T

        if (plot):
            fig = plt.figure(figsize=(16, 12))
            axs = fig.subplots(3, 3)
            labels = [r"$\alpha$", r"$\beta$", r"$\eta$"]
        
            for i, par in enumerate(samples):
                axs[2][0].set_xlabel("Ensemble step")
                axs[2][1].set_xlabel("Ensemble step")
                axs[2][2].set_xlabel("Walker Step")
                axs[i][0].set_ylabel(labels[i])
                axs[0][0].set_title("Ensemble dispersion")
                axs[0][1].set_title("Ensemble autocorrelation")
                axs[0][2].set_title("Walker mean and stdev")
                idx = np.arange(len(par))
                axs[i][0].scatter(idx, par[idx], marker='o', c='k', s=10.0, alpha=0.1, linewidth=0)
                # Get selfcorrelation using emcee
                ac = emcee.autocorr.function(par)
                
                idx = np.arange(len(ac),step=1)
                axs[i][1].scatter(idx, ac[idx], marker='o', c='k', s=10.0, alpha=0.1, linewidth=0)
                axs[i][1].axhline(alpha=1., lw=1., color='red')

                alpha_chain_mean = np.mean(alpha_chain, axis=0)
                alpha_chain_err = np.std(alpha_chain, axis=0) / np.sqrt(nwalkers)
                beta_chain_mean = np.mean(beta_chain, axis=0)
                beta_chain_err = np.std(beta_chain, axis=0) / np.sqrt(nwalkers)
                delta_chain_mean = np.mean(delta_chain, axis=0)
                delta_chain_err = np.std(delta_chain, axis=0) / np.sqrt(nwalkers)
                idx = np.arange(len(alpha_chain_mean))
                axs[0][2].errorbar(x=idx, y=alpha_chain_mean, yerr=alpha_chain_err, errorevery=50, ecolor='red',
                                   lw=0.5, elinewidth=2., color='k')
                axs[1][2].errorbar(x=idx, y=beta_chain_mean, yerr=beta_chain_err, errorevery=50, ecolor='red',
                                   lw=0.5, elinewidth=2., color='k');
                axs[2][2].errorbar(x=idx, y=delta_chain_mean, yerr=delta_chain_err, errorevery=50, ecolor='red',
                                   lw=0.5, elinewidth=2., color='k');
                
                print("Printing file:",  namemc)
                plt.savefig(namemc)
                print(namemc, "Printed")
                corner_plot(samples, labels, namecont)

        samples[:, 2] = np.exp(samples[:, 2])       
        if (nsig==1):
            a_perc = np.percentile(samples[0], [16, 50, 84]); alp =[a_perc[1], a_perc[0] - a_perc[1], a_perc[2] - a_perc[1]]
            b_perc = np.percentile(samples[1], [16, 50, 84]); bet =[b_perc[1], b_perc[0] - b_perc[1], b_perc[2] - b_perc[1]] 
            d_perc = np.percentile(samples[2], [16, 50, 84]); de =[d_perc[1], d_perc[0] - d_perc[1], d_perc[2] - d_perc[1]]
        elif(nsig==2):
            a_perc = np.percentile(samples[0], [2.3, 50, 97.7]); alp =[a_perc[1], a_perc[0] - a_perc[1], a_perc[2] - a_perc[1]]
            b_perc = np.percentile(samples[1], [2.3, 50, 97.7]); bet =[b_perc[1], b_perc[0] - b_perc[1], b_perc[2] - b_perc[1]] 
            d_perc = np.percentile(samples[2], [2.3, 50, 97.7]); de =[d_perc[1], d_perc[0] - d_perc[1], d_perc[2] - d_perc[1]]

        return [alp, bet, de]
            
    elif(gflag and (not betaflag)):
        print("")
    elif((not gflag) and bflag):
        # initial position at maximum likelihood values
        ndim = 2
        pos = [best_pars + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        # MCMC chain with 50 walkers and 1000 steps
        sampler = emcee.EnsembleSampler(nwalkers, ndim, logpost, threads=4,  args=(data,svalue,eq,gflag,bflag, moderr) )
        print("Runing MCMC ...")
        sampler.run_mcmc(pos, nsteps)
        print("Run finished")
        # Getting chains
        alpha_chain = sampler.chain[:,:,0]
        beta_chain = sampler.chain[:,:,1]
      
        # Reshaping
        alpha_chain_flat = np.reshape(alpha_chain, (nwalkers*nsteps,))
        beta_chain_flat = np.reshape(beta_chain, (nwalkers*nsteps,))
        samples = np.c_[alpha_chain_flat, beta_chain_flat].T

        if (plot):
            fig = plt.figure(figsize=(16, 8))
            axs = fig.subplots(2, 3)
            labels = [r"$\alpha$", r"$\beta$"]
            
            for i, par in enumerate(samples):
                axs[1][0].set_xlabel("Ensemble step")
                axs[1][1].set_xlabel("Ensemble step")
                axs[1][2].set_xlabel("Walker Step")
                axs[i][0].set_ylabel(labels[i])
                axs[0][0].set_title("Ensemble dispersion")
                axs[0][1].set_title("Ensemble autocorrelation")
                axs[0][2].set_title("Walker mean and stdev")
                idx = np.arange(len(par))
                axs[i][0].scatter(idx, par[idx], marker='o', c='k', s=10.0, alpha=0.1, linewidth=0)
                # Get selfcorrelation using emcee
                ac = emcee.autocorr.function(par)
                
                idx = np.arange(len(ac),step=1)
                axs[i][1].scatter(idx, ac[idx], marker='o', c='k', s=10.0, alpha=0.1, linewidth=0)
                axs[i][1].axhline(alpha=1., lw=1., color='red')
                
                alpha_chain_mean = np.mean(alpha_chain, axis=0)
                alpha_chain_err = np.std(alpha_chain, axis=0) / np.sqrt(nwalkers)
                beta_chain_mean = np.mean(beta_chain, axis=0)
                beta_chain_err = np.std(beta_chain, axis=0) / np.sqrt(nwalkers)
                idx = np.arange(len(alpha_chain_mean))
                axs[0][2].errorbar(x=idx, y=alpha_chain_mean, yerr=alpha_chain_err, errorevery=50, ecolor='red',
                                   lw=0.5, elinewidth=2., color='k')
                axs[1][2].errorbar(x=idx, y=beta_chain_mean, yerr=beta_chain_err, errorevery=50, ecolor='red',
                                   lw=0.5, elinewidth=2., color='k');
                if namemc is not None:
                    print("Printing file:",  namemc)
                    plt.savefig(namemc)
                    fig.close()
                    print(namemc, "Printed")
       
                corner_plot(samples, labels, namecont)
            
        samples[:, 2] = np.exp(samples[:, 2])
        if (nsig==1):
            a_perc = np.percentile(samples[0], [16, 50, 84]); alp =[a_perc[1], a_perc[0] - a_perc[1], a_perc[2] - a_perc[1]]
            b_perc = np.percentile(samples[1], [16, 50, 84]); bet =[b_perc[1], b_perc[0] - b_perc[1], b_perc[2] - b_perc[1]] 
        elif(nsig==2):
            a_perc = np.percentile(samples[0], [2.3, 50, 97.7]); alp =[a_perc[1], a_perc[0] - a_perc[1], a_perc[2] - a_perc[1]]
            b_perc = np.percentile(samples[1], [2.3, 50, 97.7]); bet =[b_perc[1], b_perc[0] - b_perc[1], b_perc[2] - b_perc[1]] 
        return [alp, bet]
        
    else:
        # initial position at maximum likelihood values
        ndim = 1
        pos = [best_pars + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        # MCMC chain with 50 walkers and 1000 steps
        sampler = emcee.EnsembleSampler(nwalkers, ndim, logpost, threads=4,  args=(data,svalue,eq,gflag,bflag, moderr) )
        print("Runing MCMC ...")
        sampler.run_mcmc(pos, nsteps)
        print("Run finished")
        # Getting chains
        alpha_chain = sampler.chain[:,:,0]
      
        # Reshaping
        alpha_chain_flat = np.reshape(alpha_chain, (nwalkers*nsteps,))

        samples = np.c_[alpha_chain_flat].T

        if(plot):
            fig = plt.figure(figsize=(16, 8))
            axs = fig.subplots(1, 3)
            labels = [r"$\alpha$"]
            
            for i, par in enumerate(samples):
                axs[0].set_xlabel("Ensemble step")
                axs[1].set_xlabel("Ensemble step")
                axs[2].set_xlabel("Walker Step")
                axs[0].set_ylabel(labels[i])
                axs[0].set_title("Ensemble dispersion")
                axs[1].set_title("Ensemble autocorrelation")
                axs[2].set_title("Walker mean and stdev")
                idx = np.arange(len(par))
                axs[0].scatter(idx, par[idx], marker='o', c='k', s=10.0, alpha=0.1, linewidth=0)
                # Get selfcorrelation using emcee
                ac = emcee.autocorr.function(par)
                
                idx = np.arange(len(ac),step=1)
                axs[1].scatter(idx, ac[idx], marker='o', c='k', s=10.0, alpha=0.1, linewidth=0)
                axs[1].axhline(alpha=1., lw=1., color='red')
                
                alpha_chain_mean = np.mean(alpha_chain, axis=0)
                alpha_chain_err = np.std(alpha_chain, axis=0) / np.sqrt(nwalkers)
                idx = np.arange(len(alpha_chain_mean))
                axs[2].errorbar(x=idx, y=alpha_chain_mean, yerr=alpha_chain_err, errorevery=50, ecolor='red',
                                lw=0.5, elinewidth=2., color='k')
                if namemc is not None:
                    print("Printing file:",  namemc)
                    plt.savefig(namemc)
                    fig.close()
                    print(namemc, "Printed")

                corner_plot(samples, labels, namecont)

        samples[:, 2] = np.exp(samples[:, 2])
        if (nsig==1):
            a_perc = np.percentile(samples[0], [16, 50, 84]); alp =[a_perc[1], a_perc[0] - a_perc[1], a_perc[2] - a_perc[1]]
        elif(nsig==2):
            a_perc = np.percentile(samples[0], [2.3, 50, 97.7]); alp =[a_perc[1], a_perc[0] - a_perc[1], a_perc[2] - a_perc[1]]
        return [alp]

def OneParMaxLike(best_pars,data,eq=None, svalue=None , gflag=True,bflag = True, moderr=False):    
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    from chi2 import CHI2shifted
    import numpy as np
    alpha_16 = alpha_percentil(0.16, best_pars, data, svalue=svalue,  eq=eq)
    alpha_50 = alpha_percentil(0.5, best_pars, data, svalue=svalue, eq=eq )
    alpha_84 = alpha_percentil(0.84,  best_pars, data, svalue=svalue, eq=eq)
    #print(alpha_50)
    #print(alpha_16)
    #print(alpha_84)
    print('alpha=', best_pars, '(-', alpha_50 - alpha_16, ')(+', alpha_84 - alpha_50, ')')
    
    xmin, xmax, npoints = -0.05, 0.05, 100
    alphas = np.linspace(xmin, xmax, npoints)
    #logprior and prior
    lprior = [logprior(x, gflag=gflag ,bflag = bflag)  for x in alphas]
    prior = np.exp(lprior)
    #loglikelihood and likelihood
    llike = [loglike(CHI2shifted(x,data, svalue=svalue,eq=eq, gflag=gflag, bflag=bflag, moderr=moderr ))  for x in alphas]
    like= np.exp(llike)
    #logposterior and posterior
    lpost = [logpost(x, data,svalue=svalue, eq=eq, gflag=gflag, bflag=bflag, moderr=moderr ) for x in alphas]
    posterior = np.exp(lpost)
    plt.clf()
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"Prior")
    plt.xlim([xmin, xmax])
    plt.plot(alphas, prior)
    print("Printing file: Prior_only_alpha.png ")
    plt.savefig('Prior_only_alpha.png')#, dpi=150)
    
    plt.clf()
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"Likelihood")
    plt.xlim([xmin, xmax])
    plt.plot(alphas, like)
    print("Printing file: Likelihood_only_alpha.png ")
    plt.savefig('Likelihood_only_alpha.png')#, dpi=150)
    
    #plt.clf()
    #plt.xlabel(r"$\alpha$")
    #plt.ylabel(r"Posterior probability")
    #plt.xlim([xmin, xmax])
    #print(lpost)
    #plt.plot(alphas, posterior)
    #plt.axvline(x=alpha_16, linestyle="dashed", color="red")
    #plt.axvline(x=alpha_84, linestyle="dashed", color="red")
    #plt.axvline(x=alpha_50, color="red");
    #print("Printing file: Posterior_only_alpha.png ")
    #plt.savefig('Posterior_only_alpha.png')#, dpi=150)

        
        




   
        
