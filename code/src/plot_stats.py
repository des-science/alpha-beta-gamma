import matplotlib.pyplot as plt

#plt.style.use('supermongo.mplstyle')
plt.style.use('SVA1StyleSheet.mplstyle')

def pretty_rho0(meanr, rho, sig, sqrtn):
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0]/sqrtn, color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0]/sqrtn, color='blue', ls='', marker='o')
    rho0_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    plt.legend([rho0_line],[r'$\rho_0(\theta)$'],loc='upper right', fontsize=24)
    #plt.ylim( [1.e-9, 5.e-6] )
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\rho(\theta)$', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()

def pretty_rho1(meanr, rho, sig, sqrtn, rho3=None, sig3=None, rho4=None, sig4=None):
    import matplotlib.patches as mp
    
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0]/sqrtn, color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0]/sqrtn, color='blue', ls='', marker='o')
    rho1_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    if rho3 is not None:
        plt.plot(meanr*1.03, rho3, color='green')
        plt.plot(meanr*1.03, -rho3, color='green', ls=':')
        plt.errorbar(meanr[rho3>0]*1.03, rho3[rho3>0], yerr=sig3[rho3>0]/sqrtn, color='green', ls='', marker='s')
        plt.errorbar(meanr[rho3<0]*1.03, -rho3[rho3<0], yerr=sig3[rho3<0]/sqrtn, color='green', ls='', marker='s')
        rho3_line = plt.errorbar(-meanr, rho3, yerr=sig3, color='green', marker='s')
    if rho4 is not None:
        plt.plot(meanr*1.06, rho4, color='red')
        plt.plot(meanr*1.06, -rho4, color='red', ls=':')
        plt.errorbar(meanr[rho4>0]*1.06, rho4[rho4>0], yerr=sig4[rho4>0]/sqrtn, color='red', ls='', marker='^')
        plt.errorbar(meanr[rho4<0]*1.06, -rho4[rho4<0], yerr=sig4[rho4<0]/sqrtn, color='red', ls='', marker='^')
        rho4_line = plt.errorbar(-meanr, rho4, yerr=sig4, color='red', marker='^')
    #sv_req = mp.Patch(color='#FFFF82')
    if rho3 is not None and rho4 is not None:
        plt.legend([rho1_line, rho3_line, rho4_line],
                   [r'$\rho_1(\theta)$', r'$\rho_3(\theta)$', r'$\rho_4(\theta)$'],
                   loc='upper right', fontsize=24)
        #plt.ylim( [1.e-9, 5.e-6] )
        #plt.ylim( [1.e-9, 2.e-5] )
        plt.ylim( [5.e-9, 5.e-6] )
    elif True:
        plt.legend([rho1_line, sv_req],
                   [r'$\rho_1(\theta)$', r'Requirement'],
                   loc='upper right')
        plt.ylim( [1.e-9, 5.e-6] )
    else: # For talk
        plt.legend([rho1_line, sv_req],
                   [r'$\rho_1(\theta)$',
                    r'Requirements for $d\sigma_8/\sigma_8 < 0.03$'],
                   loc='upper right')
        plt.ylim( [1.e-9, 3.e-6] )
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\rho(\theta)$', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()

def pretty_rho2(meanr, rho, sig, sqrtn, rho5=None, sig5=None ):
    import matplotlib.patches as mp
    # The requirements on rho2 are less stringent.  They are larger by a factor 1/alpha.
    # Let's use alpha = 0.03.
    alpha = 0.03

    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0]/sqrtn, color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0]/sqrtn, color='blue', ls='', marker='o')
    rho2_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    if rho5 is not None:
        plt.plot(meanr*1.03, rho5, color='green')
        plt.plot(meanr*1.03, -rho5, color='green', ls=':')
        plt.errorbar(meanr[rho5>0]*1.03, rho5[rho5>0], yerr=sig5[rho5>0]/sqrtn, color='green', ls='', marker='s')
        plt.errorbar(meanr[rho5<0]*1.03, -rho5[rho5<0], yerr=sig5[rho5<0]/sqrtn, color='green', ls='', marker='s')
        rho5_line = plt.errorbar(-meanr, rho5, yerr=sig5, color='green', marker='s')
    #sv_req = mp.Patch(color='#FFFF82')

    if rho5 is not None :
        plt.legend([rho2_line, rho5_line ],
                   [r'$\rho_2(\theta)$', r'$\rho_5(\theta)$'],
                   loc='upper right', fontsize=24)
        #plt.ylim( [1.e-7, 5.e-4] )
        plt.ylim( [5.e-8, 1.e-5] )
    elif True: # For paper
        plt.legend([rho2_line, sv_req],
                   [r'$\rho_2(\theta)$', r'Requirement'],
                   loc='upper right')
        plt.ylim( [1.e-7, 5.e-4] )
    else:
        plt.legend([rho2_line, sv_req],
                   [r'$\rho_2(\theta)$',
                    r'Requirements for $d\sigma_8/\sigma_8 < 0.03$'],
                   loc='upper right')
        plt.ylim( [1.e-7, 3.e-4] )
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(r'$\rho(\theta)$', fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()

def pretty_sigma(meanr, rho, sig, sqrtn, title):
    plt.plot(meanr, rho, color='blue')
    plt.plot(meanr, -rho, color='blue', ls=':')
    plt.errorbar(meanr[rho>0], rho[rho>0], yerr=sig[rho>0]/sqrtn, color='blue', ls='', marker='o')
    plt.errorbar(meanr[rho<0], -rho[rho<0], yerr=sig[rho<0]/sqrtn, color='blue', ls='', marker='o')
    rho0_line = plt.errorbar(-meanr, rho, yerr=sig, color='blue', marker='o')
    plt.legend([rho0_line],[title],loc='upper right', fontsize=24)
    #plt.ylim( [1.e-9, 5.e-6] )
    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.xlim( [0.5,300.] )
    plt.xlabel(r'$\theta$ (arcmin)', fontsize=24)
    plt.ylabel(title, fontsize=24)
    plt.xscale('log')
    plt.yscale('log', nonposy='clip')
    plt.tight_layout()
