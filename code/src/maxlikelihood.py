#Log natural of the prior 
def logprior(pars, gflag=True,bflag = True):
    import numpy as np
    if(gflag and betaflag):
        #print("Using alpha, beta and gamma")
        alpha, beta, gamma = pars
    elif(gflag and (not betaflag)):
        #print("Using alpha and gamma")
        alpha, gamma = pars
    elif((not gflag) and bflag):
        #print("Using alpha and beta")
        alpha, beta = pars
        if (-3.<alpha<3)and(-3< beta<5):
            return 0.0
        return -np.inf
    else:
        #print("Using only alpha")
        alpha = pars
        if ( - 0.1 < alpha < 0.1):
            return 0.0
        return -np.inf
##Log natural of the likelihood function. Gaussian.
def loglike(chisq):
    return -0.5*chisq
##Log natural of the posterior
def logpost(pars, data, eq=None, gflag=True,bflag = True):
    import numpy as np
    from chi2 import CHI2
    chisq = CHI2(pars, data,eq=eq, gflag=gflag, bflag=bflag )
    lp = logprior(pars, gflag=gflag,bflag = bflag)
    if not np.isfinite(lp):
        return -np.inf
    return lp + loglike(chisq)

def MCMC(best_pars,data, eq=None, gflag=True,bflag = True):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    import emcee
    import corner
    import numpy as np
    if(gflag and betaflag):
        print("")
    elif(gflag and (not betaflag)):
        print("")
    elif((not gflag) and bflag):
        # initial position at maximum likelihood values
        ndim, nwalkers, nsteps = 2, 50, 1000
        pos = [best_pars + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        # MCMC chain with 50 walkers and 1000 steps
        sampler = emcee.EnsembleSampler(nwalkers, ndim, logpost, threads=4,  args=(data,eq,gflag,bflag) )
        print("Runing MCMC ...")
        sampler.run_mcmc(pos, nsteps)
        print("Run finished")
        # Getting chains
        alpha_chain = sampler.chain[:,:,0]
        beta_chain = sampler.chain[:,:,1]
        # Average and standard deviation between chains
        alpha_chain_mean = np.mean(alpha_chain, axis=0)
        alpha_chain_std = np.std(alpha_chain, axis=0) / np.sqrt(nwalkers)
        # Reshaping
        alpha_chain_flat = np.reshape(alpha_chain, (nwalkers*nsteps,))
        beta_chain_flat = np.reshape(beta_chain, (nwalkers*nsteps,))

        fig = plt.figure(figsize=(16, 8))
        axs = fig.subplots(2, 3)
        labels = [r"$\alpha$", r"$\beta$"]
        samples = np.c_[alpha_chain_flat, beta_chain_flat].T
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
        print("Printing file: mcmc_alpha_beta.pdf")
        plt.savefig('mcmc_alpha_beta.pdf')
        print("mcmc_alpha_beta.pdf Printed")

        # corner plot
        burn = 5000
        samples_burned = np.c_[[par[burn:] for par in samples]]
        fig = corner.corner(samples_burned.T, labels=labels,
                    quantiles=[0.16, 0.5, 0.84], 
                    levels=(1-np.exp(-0.5), 1-np.exp(-2), 1-np.exp(-9./2)), #1sigma, 2sigma and 3sigma contours
                    show_titles=True, title_kwargs={"fontsize": 12},
                    smooth1d=None, plot_contours=True,
                    no_fill_contours=False, plot_density=True,)
        print("Printing file: contour_alpha_beta.pdf")
        plt.savefig('contour_alpha_beta.pdf')
        print("contour_alpha_beta.pdf Printed")


    else:
        print("")

        
    
    
def alpha_percentil(p, mu, data, eq=None):
    from scipy.integrate import quad
    from scipy.optimize import brentq
    import numpy as np
    max_loglike= logpost(mu,data,eq=None, gflag=False, bflag=False)
    print(mu,  max_loglike)
    def cumulative_like(alpha):
        integral, _ = quad(lambda x: np.exp(logpost(x, data, eq=eq, gflag=False, bflag=False ) - max_loglike), -np.inf, alpha)
        return integral
    total_int = cumulative_like(np.inf)
    #print(total_int)
    perc = brentq(lambda x: cumulative_like(x) - p*total_int, -np.inf, np.inf)
    return perc

def all_posterior_info(best_pars,data, eq=None,  gflag=True, bflag = True):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('SVA1StyleSheet.mplstyle')
    import numpy as np
    from chi2 import CHI2
    if(gflag and betaflag):
        #print("Using alpha, beta and gamma")
        alpha, beta, gamma = pars
    elif(gflag and (not bflag)):
        #print("Using alpha and gamma")
        alpha, gamma = pars
    elif((not gflag) and bflag):
        #print("Using alpha and beta")
        MCMC(best_pars,data, eq=eq, gflag=gflag,bflag = bflag)
    else:
        #print("Using only alpha")
        alpha_16 = alpha_percentil(0.16, best_pars, data, eq=eq)
        alpha_50 = alpha_percentil(0.5, best_pars, data, eq=eq)
        alpha_84 = alpha_percentil(0.84,  best_pars, data, eq=eq)
        #print(alpha_50)
        #print(alpha_16)
        #print(alpha_84)
        print('alpha=', best_pars, '(-', alpha_50 - alpha_16, ')(+', alpha_84 - alpha_50, ')')
        
        
        xmin, xmax, npoints = -0.5, 0.5, 100
        alphas = np.linspace(xmin, xmax, npoints)
        #logprior and prior
        lprior = [logprior(x, gflag=gflag ,bflag = bflag)  for x in alphas]
        prior = np.exp(lprior)
        #loglikelihood and likelihood
        llike = [loglike(CHI2(x, data,eq=eq, gflag=gflag, bflag=bflag ))  for x in alphas]
        like= np.exp(llike)
        #logposterior and posterior
        lpost = [logpost(x, data, eq=eq, gflag=gflag, bflag=bflag ) for x in alphas]
        posterior = np.exp(lpost)
        plt.clf()
        plt.xlabel(r"$\alpha$")
        plt.ylabel(r"Prior")
        plt.xlim([xmin, xmax])
        plt.plot(alphas, prior)
        print("Printing file: Prior_only_alpha.pdf ")
        plt.savefig('Prior_only_alpha.pdf')#, dpi=150)

        plt.clf()
        plt.xlabel(r"$\alpha$")
        plt.ylabel(r"Likelihood")
        plt.xlim([xmin, xmax])
        #plt.plot(alphas, posterior)
        #plt.axvline(x=alpha_16, linestyle="dashed", color="red")
        #plt.axvline(x=alpha_84, linestyle="dashed", color="red")
        #plt.axvline(x=alpha_50, color="red");
        #print("Printing file: Posterior_only_alpha.pdf ")
        #plt.savefig('Posterior_only_alpha.pdf')#, dpi=150)

        #plt.clf()
        #plt.xlabel(r"$\alpha$")
        #plt.ylabel(r"Posterior probability")
        #plt.xlim([xmin, xmax])
        #print(lpost)
        #plt.plot(alphas, posterior)
        #plt.axvline(x=alpha_16, linestyle="dashed", color="red")
        #plt.axvline(x=alpha_84, linestyle="dashed", color="red")
        #plt.axvline(x=alpha_50, color="red");
        #print("Printing file: Posterior_only_alpha.pdf ")
        #plt.savefig('Posterior_only_alpha.pdf')#, dpi=150)




   
        
