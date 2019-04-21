def read_corr(stat_file):
    import numpy as np
    import fitsio
    alldata =  fitsio.read(stat_file, ext=2)
    covmat =  fitsio.read(stat_file, ext=1)
    angdist = alldata['ANG']; corr =  alldata['VALUE']
    return  angdist, corr,  covmat

