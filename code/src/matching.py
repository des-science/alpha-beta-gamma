import os
import pandas 

def parse_args():
    import argparse
    from astropy import units as u
    parser = argparse.ArgumentParser(description='A program to match two files using their RA a DEC ')
    
    parser.add_argument('--file', default='',
                        help='File to be modified, the output will be this file plus other additional field')
    parser.add_argument('--fileref', default='',
                        help='File of reference, this file have information to be added to the analysis file ')
    parser.add_argument('--filename', default='out.fits',
                        help='Name for the output file ')
    parser.add_argument('--d2d', default='', type=float, 
                        help='Maximum angular separation to be accept as a rigth match (in units of degree) ')
    #parser.add_argument('--unit', default= , type=u, 
    #                    help='RA and DEC units in the files and for the angular separation (u.degre) ')
    parser.add_argument('--fields', nargs='+',  type=str, 
                        help='Fields from fileref to be include in file, --fields string1 string2 string3 ... ')
    parser.add_argument('--overwrite', default=False, action='store_const', const=True,
                        help='boolean to overwrite the field in the original file')
    parser.add_argument('--format', default='fits', type=str, 
                        help='format of the input files, by the moment txt or fits, instert it as string')

    args = parser.parse_args()
    return args

def checknames(fileaux):
    import fitsio
    from fitsio import FITS
    fit = fitsio.FITS(fileaux)
    candidate =  ['Ra', 'RA', 'ra', 'Dec', 'DEC',  'dec']
    names = []
    for name in candidate:
        if name in fit[1].get_colnames():
            names.append(name)
    if (len(names) != 2) :
        print(' The file does not have RA and DEC or their names are not in the list:',  candidate)
    return names
        
def testinputfiles(args):
    if args.file != '':
        print ('\n File ',args.file, 'sucessfully read. OK! \n' )
    else:
        print ('No input file')

    if args.fileref != '':
        print ('\n Fileref ',args.fileref, 'sucessfully read. OK! \n' )
    else:
        print ('No input file')

def match_objects(df1, df2,  names,  args):
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    import numpy as np
    import pandas as p

    c = SkyCoord(ra=list(df1[names[0]])* u.degree, dec=list(df1[names[1]])*u.degree)  
    catalog = SkyCoord(ra=list(df2[names[2]])*u.degree, dec=list(df2[names[3]])*u.degree)
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)  
    matches = catalog[idx]

    df2 =  df2.reindex(index=idx)
    #now df2 starts with the rows having the information of the match
    # and nothing about the other rows of df2
    if (len(matches)== len(df1)):
        df1['d2d'] = list(d2d.degree)
        #df1['cmt'] =  list(df2.loc[df2.index.isin(idx) , 'CM_T' ])
        newcols = {}
        for name in args.fields:
            newcols[name] =  list(df2.loc[df2.index.isin(idx) , name ])
            df1[name] = p.Series(newcols[name],  index = df1.index )
            boolean =  ( df1['d2d'] >  args.d2d )
            df1.loc[boolean, name] = None
    return df1
       
def match_objects_test1(df1, df2):
    import numpy as np
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    dat1 = np.loadtxt(df1,  dtype={'names': ('ra', 'dec'), 'formats': ('f4', 'f4')})
    dat2 =  np.loadtxt(df2,  dtype={'names': ('RA', 'DEC', 'EXT_MOF'), 'formats': ('f4', 'f4',  'f4')})
    print('\n Data 1: \n')
    PrettyPrint(dat1)
    print('\n Data 2: \n')
    PrettyPrint(dat2)
    
    c = SkyCoord(ra=list(dat1['ra'])* u.degree, dec=list(dat1['dec'])*u.degree)  
    catalog = SkyCoord(ra=list(dat2['RA'])*u.degree, dec=list(dat2['DEC'])*u.degree)
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)  
    matches = catalog[idx]
 
    print('\n IDX:', idx)
    #print(d2d.degree)
    
    names = dat1.dtype.names
    dt = dat1.dtype
    dt = dt.descr + [('ext_mof', '<f8'),  ('d2d', '<f8')]
    out_data =  np.recarray(shape= (len(dat1),) , dtype=dt)
    for name in names:
        out_data[name] = dat1[name]

    if (len(matches)== len(dat1)):
        for i in range (0, len(matches)):
            out_data['d2d'][i] = d2d.degree[i]
            out_data['ext_mof'][i] =  dat2['EXT_MOF'][idx[i]]
    print('\n Final data 1: \n')
    PrettyPrint(out_data)    
    #np.savetxt('out_data.txt', out_data, fmt='%f')
def match_objects_test2(arg1, arg2):
    import numpy as np
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    import fitsio
    dat1 = np.loadtxt(arg1,  dtype={'names': ('ra', 'dec'), 'formats': ('f4', 'f4')})
    dat1 = dat1.astype(dat1.dtype.newbyteorder('='))
    df1 = pandas.DataFrame(dat1)
    dat2 =  np.loadtxt(arg2,  dtype={'names': ('RA', 'DEC', 'EXT_MOF'), 'formats': ('f4', 'f4',  'f4')})
    dat2 = dat2.astype(dat2.dtype.newbyteorder('='))
    df2 = pandas.DataFrame(dat2)
    print('\n Data 1: \n', df1)
    print('\n Data 2: \n', df2)

    
    c = SkyCoord(ra=list(dat1['ra'])* u.degree, dec=list(dat1['dec'])*u.degree)  
    catalog = SkyCoord(ra=list(dat2['RA'])*u.degree, dec=list(dat2['DEC'])*u.degree)
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)  
    matches = catalog[idx]
 
    print('\n IDX:',  idx)
    #print(d2d.degree)

    df2 =  df2.reindex(index=idx)

    print('\n Data 2 reindexed: \n', df2)
    if (len(matches)== len(dat1)):
        df1['d2d'] = list(d2d.degree)
        df1['ext_mof'] =  list(df2.loc[df2.index.isin(idx) ,  'EXT_MOF'])

    print('\n Final data 1: \n', df1)
    
    '''
    for i in range (0, len(matches)):
        out_data['d2d'][i] = d2d.degree[i]
        out_data['ext_mof'][i] =  dat2['EXT_MOF'][idx[i]]
    np.savetxt('out_data.txt', df1, fmt='%f')
    fitsio.write('outdatas.fits', df1.to_records(index=False), clobber=True)
    '''

def PrettyPrint(data):
    from prettytable import PrettyTable
    x = PrettyTable(data.dtype.names)    
    for row in data:
        x.add_row(row)
        # Change some column alignments; default was 'c'
        #x.align['column_one'] = 'r'
        #x.align['column_3'] = 'l'
    print(x)
    
def write_fit(data, file_name):
    import fitsio
    print("Writing data to ",file_name)
    fitsio.write(file_name, data.to_records(index=False), clobber=True)
    #fitsio.write(file_name, data, clobber=True)
    
def main():
    '''
    from astropy.io import fits
    hdu1=  fits.open(args.file)
    hdu2 =  fits.open(args.fileref)
    data1 = hdu1[1].data
    data2 = hdu2[1].data
    '''
   
    import fitsio
    import pandas

    args = parse_args()
    testinputfiles(args)

    data1 =  fitsio.read(args.file)
    data1 = data1.astype(data1.dtype.newbyteorder('='))
    df1 = pandas.DataFrame(data1)

    data2 =  fitsio.read(args.fileref)
    data2 = data2.astype(data2.dtype.newbyteorder('='))
    df2 = pandas.DataFrame(data2)

    names1 =  checknames(args.file)
    names2 =  checknames(args.fileref)
    names = names1 + names2 
 
    outdata = match_objects(df1, df2,  names,  args)
    write_fit(outdata, args.filename)

    

if __name__ == "__main__":
    main()
    




    
    
 
