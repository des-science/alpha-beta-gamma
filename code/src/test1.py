def main():
    from readjson import read_rhos
    usualrhos =  "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/used/rho_all_reserved_irz_usual.json"
    modrhos = "/home2/dfa/sobreira/alsina/catalogs/output/alpha-beta-gamma/used/rho_all_reserved_irz_modified.json"
    umeanr, urho0p, urho1p, urho2p, urho3p, urho4p, urho5p, usig_rho0, usig_rho1, usig_rho2, usig_rho3, usig_rho4, usig_rho5 = read_rhos(usualrhos)
    meanr, rho0p, rho1p, rho2p, rho3p, rho4p, rho5p, sig_rho0, sig_rho1, sig_rho2, sig_rho3, sig_rho4, sig_rho5 = read_rhos(modrhos)

    print(rho1p - urho1p)
    
if __name__ == "__main__":
    main()
