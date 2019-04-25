def read_corr(stat_file, maxscale=None,  maxbin=None):
    import numpy as np
    import fitsio
    alldata =  fitsio.read(stat_file, ext=2)
    covmat =  fitsio.read(stat_file, ext=1)
    angdist = alldata['ANG']; corr =  alldata['VALUE']
    if maxbin is not None:
        angdist = angdist[:maxbin]
        corr = corr[:maxbin]
        covmat = covmat[:maxbin,:maxbin]
    if maxscale is not None:
        angdist = angdist[angdist<maxscale]
        idx = len(angdist)
        corr= corr[:idx]
        covmat = covmat[:idx,:idx]
    return  angdist, corr,  covmat

