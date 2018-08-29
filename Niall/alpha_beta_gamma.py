import yaml
import corrtools
import collections
import twopoint
import argparse
from os.path import join as pj
from numpy.lib.recfunctions import append_fields
import numpy as np
import treecorr
import fitsio
import os
import sys
import pylab
#import mpi4py
#from mpi4py import MPI
#from mpi4py.MPI import ANY_SOURCE
from corrtools import sample_cov
from corrtools.y1_tools import *
import pickle
import nmbot
import healpix_util as hu
import healpy as hp
#homogenised names
homog_colnames={'e1':'e1', 'e2':'e2', 'ra':'ra', 'dec':'dec', 'psf_e1':'psf_e1', 'psf_e2':'psf_e2','weight':'weight'}
shape_cols=['e1','e2','psf_e1','psf_e2','flags','weight','snr','size','psf_size']

def read_psf_data(args):
    columns = ['ra', 'dec', 'e1', 'e2', 'de1', 'de2', args.psf_model_prefix+'_e1', args.psf_model_prefix+'_e2']
    #columns = ['ra', 'dec', args.psf_model_prefix+'_e1', args.psf_model_prefix+'_e2','de1', 'de2']
    if args.star_weight_col:
        columns.append(args.star_weight_col)
    if args.dTp or args.dTpT:
        columns+=['dsize',args.psf_model_prefix+'_size']
    psf_data=fitsio.read(args.psf_cat, columns=columns)
    #good_psf=(psf_data['psfex_e1']>-1.)*(psf_data['psfex_e1']<1.)*(psf_data['e1']>-1.)*(psf_data['e1']<1.)
    #print 'keeping %d/%d psfs'%(good_psf.sum(),len(psf_data))
    #psf_data=psf_data[good_psf]
    if args.star_cut_dict:
        print 'applying cut file %s to star catalog'%args.star_cut_dict
        with open(args.star_cut_dict, 'rb') as f:
            cut_dict=yaml.load(f)
        psf_data,_ = nmbot.apply_cut_dict(psf_data, cut_dict, verbose=True)
    return psf_data

def zero_fill(zeros):
    return {'xip':zeros, 'xim':zeros, 'theta':zeros}

def parse_args():

    description = ""
    parser = argparse.ArgumentParser(description=description, add_help=True,
                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('shape_cat', type=str)
    parser.add_argument('psf_cat',type=str)
    parser.add_argument('outdir',type=str)
    parser.add_argument('--cut_dict', type=str)
    parser.add_argument('--z_bin_lims', type=float, nargs='*', default=None)
    parser.add_argument('--star_weight_col',type=str,default=None)
    parser.add_argument('--pipeline',type=str,default='metacal')
    parser.add_argument('--bin_slop',type=float,default=0.5)
    parser.add_argument('--dTp',action='store_true',default=False)
    parser.add_argument('--dTpT',action='store_true',default=False)
    parser.add_argument('--psf_model_prefix',default='psf',help="prefix to psf-model columns e.g. should be 'psf' if psf-model e1 columns is called 'psf_e1'")
    parser.add_argument('--star_cut_dict',default=None, help='cut dict for star catalog')
    parser.add_argument('--gold_fp_mask',default=None)
    parser.add_argument('--gold_br_mask',default=None)
    parser.add_argument('--nbins',default=15,type=int)
    parser.add_argument('--q_cat',type=str, default=None)
    parser.add_argument('--q_coeff_cat',type=str, default=None)
    parser.add_argument('--test_nrows',type=int,default=None)
    return parser.parse_args()

def main():

    args=parse_args()
    if args.q_cat is not None:
        assert args.q_coeff_cat is not None

    psf_data=read_psf_data(args)

    #if requested, apply gold footprint mask and/or gold badregions mask
    use = np.ones(len(psf_data), dtype=bool)
    if args.gold_fp_mask:
        fp=hu.Map("ring",hp.read_map(args.gold_fp_mask))
        fp_vals=fp.get_mapval(psf_data['ra'],psf_data['dec'])
        use[fp_vals<1]=False
        print float(use.sum())/len(use)
    if args.gold_br_mask:
        br=hu.Map("ring",hp.read_map(args.gold_br_mask))
        br_vals=br.get_mapval(psf_data['ra'],psf_data['dec'])
        use[br_vals>0]=False
        print float(use.sum())/len(use)
    if use.sum()!=len(use):
        print 'gold masks leave fraction %f of stars'%(float(use.sum())/len(use))
    psf_data=psf_data[use]

    theta_min,theta_max,nbins=0.25,250.,args.nbins
    bin_slop=args.bin_slop
    sep_units='arcmin'
    gg=treecorr.GGCorrelation(nbins=nbins, min_sep=theta_min, max_sep=theta_max, sep_units=sep_units,verbose=1,bin_slop=bin_slop)

    if args.z_bin_lims:
        num_z_bins = len(args.z_bin_lims)-1
    else:
        num_z_bins = 1

    shape_colnames=homog_colnames

    #read data
    if args.pipeline=='metacal':
        shape_data, shape_mask_0, shape_masks_sel, common_mask, read_rows_union, area_deg = get_mcal_cat(args.shape_cat, 
                                                                                               args.cut_dict, 
                                                                                               mcal_cols=['e1','e2','psf_e1','psf_e2','R11','R22','ra','dec'], 
                                                                                               test_nrows=args.test_nrows,
                                                                                               add_Rmean_col=True)

        if args.z_bin_lims:
            #supplying z_bin_lims implies you want to some 'tomographic' redshift bin splitting, so read in some redshift info too
            z_arrays = get_mcal_photoz(read_rows=read_rows_union)
            #pylab.hist(z_arrays[0],bins=20)
            #pylab.show()
            #first one is the one used for the unsheared binning (and therefore the measurement), so add this to the shape_data array
            shape_data = nmbot.add_field(shape_data, [("mean_z",float)], [z_arrays[0]])

        #now make a McalCat and apply correction
        shape_cat = corrtools.McalCat(shape_data, shape_mask_0, shape_masks_sel, quantities=[('e1','e2'),('psf_e1','psf_e2')])
        print len(shape_cat.arr_data)
        if args.z_bin_lims:
            print shape_cat.arr_data['mean_z'].min(), shape_cat.arr_data['mean_z'].max()
            shape_cat.redshift_split(zcol='mean_z', z_bin_lims=args.z_bin_lims)
            shape_cat.sheared_redshift_split(z_arrs_sheared = z_arrays[1:])
            shape_cat.apply_mcal_correction()
            print shape_cat.bin_stats

    else:
        assert args.pipeline=="im3shape"
        shape_data, read_rows, area_deg = get_im3_cat(args.shape_cat, args.cut_dict, im3_cols=IM3_COLS_DEFAULT+['psf_e1','psf_e2'],
                                            test_nrows=None, apply_c=True)
        print shape_data.dtype.names
        if args.z_bin_lims:
            z_array = get_photoz(read_rows=read_rows)
            shape_data = nmbot.add_field(shape_data, [("mean_z",float)], [z_array])
        shape_cat = corrtools.GalCat(shape_data, quantities=[('e1','e2'),('psf_e1','psf_e2')])
        if args.z_bin_lims:
            shape_cat.redshift_split(zcol='mean_z', z_bin_lims=args.z_bin_lims)
            shape_cat.apply_m_correction()

        print shape_cat.bin_stats        
        weight_col='weight'

    if args.q_cat:
        q_data = fitsio.read(args.q_cat, rows=read_rows_union)[shape_mask_0]
        if args.z_bin_lims:
            q_data = q_data[(shape_data[shape_mask_0]["mean_z"]>args.z_bin_lims[0])*(shape_data["mean_z"][shape_mask_0]<args.z_bin_lims[-1])]
        print "len(q_data), len(shape_cat.arr_data)",len(q_data), len(shape_cat.arr_data)
        #shape_data,q_data=shape_data[shape_mask_0],q_data[shape_mask_0]
        #shape_data = nmbot.add_field(shape_data, [('de1',float),('de2',float)], [q_data['de1'],q_data['de2']])
        use = np.ones(len(q_data),dtype='bool')
        use[np.isnan(q_data['e1'])]=False
        use[q_data['de1']<-1.]=False
        use[q_data['de2']<-1.]=False
        print 'fraction %f has bad q data, will use mean qs for these objets'%(1 - float(use.sum())/len(use))
        mean_q1 = q_data['de1'][use].mean()
        mean_q2 = q_data['de2'][use].mean()
        q_data['de1'][~use] = mean_q1
        q_data['de2'][~use] = mean_q2
        #now compute corrections
        with open(args.q_coeff_cat,'r') as f:
            coeff_lines=f.readlines()
        for l in coeff_lines:
            if l[0]=='#':
                continue
            l_entries = (l.strip()).split()
            zbin, x, y, alpha = int(l_entries[0]), l_entries[1], l_entries[2], float(l_entries[3])
            if (x=='de1' and y=='e1'):
                shape_cat.arr_data['e1'][shape_cat.zbin_masks[zbin]] -= alpha * q_data['de1'][shape_cat.zbin_masks[zbin]]
            elif (x=='de2' and y=='e2'):
                shape_cat.arr_data['e2'][shape_cat.zbin_masks[zbin]] -= alpha * q_data['de2'][shape_cat.zbin_masks[zbin]]
            else:
                continue
            print 'zbin, x, y, alpha',zbin, x, y, alpha

    #sh_cat,im3_psf_cat = make_shape_treecorr_cats(shape_data, shape_colnames)
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    #recalculate post-correction bin stats and save
    shape_cat.get_bin_stats()
    with open(pj(args.outdir, "bin_stats_shape.pkl"),"wb") as f:
        pickle.dump(shape_cat.bin_stats, f)
              
    print psf_data.dtype.names
    psf_corr_names = ['P','p','q']
    #psf_quantities = [(args.psf_model_prefix+'_e1', args.psf_model_prefix+'_e2'),('de1','de2')]
    psf_quantities = [('e1', 'e2'),(args.psf_model_prefix+'_e1', args.psf_model_prefix+'_e2'),('de1','de2')]
    if args.dTpT:
        psf_data=nmbot.add_field(psf_data, [('dse_1',float),('dse_2',float)], 
                                 [psf_data['dsize']*psf_data['e1']/psf_data['size'], 
                                  psf_data['dsize']*psf_data['e2']/psf_data['size']])
        psf_quantities.append(('dse_1','dse_2'))
        psf_corr_names.append('dse')

    elif args.dTp:
        psf_data=nmbot.add_field(psf_data, [('dse_1',float),('dse_2',float)], 
                                 [psf_data['dsize']*psf_data['e1'], 
                                  psf_data['dsize']*psf_data['e2']])
        psf_quantities.append(('dse_1','dse_2'))
        psf_corr_names.append('dse')

                                 
    psf_cat = corrtools.GalCat(psf_data, x_col='ra', y_col='dec', quantities=psf_quantities, w_col=args.star_weight_col)
    psf_cat.redshift_split()
    with open(pj(args.outdir, "bin_stats_psf.pkl"),"wb") as f:
        pickle.dump(psf_cat.bin_stats, f)
    
    sh = (shape_colnames['e1'],shape_colnames['e2'])
    sh_psf = (shape_colnames['psf_e1'],shape_colnames['psf_e2'])
    sh_quantities = [sh, sh_psf]

    cross_corrs=[]
    cross_specs=[]
    psf_auto_specs=[]
    varxi_arr = []
    #first do shear cross psf(gal)
    epg = corrtools.CorrMeasurement(shape_cat, sh, q2=sh_psf, XX=gg, sub_mean=True)
    spec_epg,varxis = epg.get_spec(['epg']*2, 'NZ_DUMMY', kernel_name2='NZ_DUMMY', ret_varxi_arrs=True)
    cross_corrs.append(epg)
    cross_specs.append(spec_epg[0])
    varxi_arr+=list(varxis[0].copy())
    #and psf(gal) auto
    #do gp auto
    pg = corrtools.CorrMeasurement(shape_cat, sh_psf, XX=gg, sub_mean=True)
    spec_pg, varxis = pg.get_spec(['pgpg']*2, 'NZ_DUMMY', ret_varxi_arrs=True)
    varxi_arr+=list(varxis[0].copy())
    psf_auto_specs.append(spec_pg[0])

    for i in range(len(psf_quantities)):
        q2_i, name_i = psf_quantities[i], psf_corr_names[i]
        # e x p
        ep = corrtools.CorrMeasurement(shape_cat, sh, gal_cat2=psf_cat, q2=q2_i, XX=gg, sub_mean=True)
        sp, varxis = ep.get_spec(['e'+name_i]*2, 'NZ_DUMMY', kernel_name2='NZ_DUMMY', ret_varxi_arrs=True)
        cross_corrs.append(ep)
        cross_specs.append(sp[0]) #just use xip
        varxi_arr += list(varxis[0].copy())
        #do p(gal) x p
        pgp = corrtools.CorrMeasurement(shape_cat, sh_psf, gal_cat2=psf_cat, q2=q2_i, XX=gg, sub_mean=True)
        spec_pgp, varxis=pgp.get_spec(['pg'+name_i]*2, 'NZ_DUMMY', kernel_name2='NZ_DUMMY', ret_varxi_arrs=True)
        psf_auto_specs.append(spec_pgp[0])
        varxi_arr += list(varxis[0].copy())
        #sp = corrtools.SpectrumMeasurement.from_galcats(['e'+name_i]*2, shape_cat, sh, 'NZ_DUMMY', gal_cat2=psf_cat, q2=q2_i, 
        #                                                XX=gg, kernel_name2='NZ_DUMMY')
        for j in range(i,len(psf_quantities)):
            print name_i,psf_corr_names[j]
            if i==j:
                pp = corrtools.CorrMeasurement(psf_cat, q2_i, XX=gg, sub_mean=True)
                sp,varxis = pp.get_spec([name_i+name_i]*2, 'NZ_DUMMY', ret_varxi_arrs=True)
                #pp = corrtools.SpectrumMeasurement.from_galcats([name_i+name_i]*2, psf_cat, q2_i, 'NZ_DUMMY', XX=gg)
            else:
                q2_j, name_j = psf_quantities[j], psf_corr_names[j]
                pp = corrtools.CorrMeasurement(psf_cat, q2_i, XX=gg, q2=q2_j, sub_mean=True)
                sp,varxis = pp.get_spec([name_i+name_j]*2, 'NZ_DUMMY', ret_varxi_arrs=True)
                #pp = corrtools.SpectrumMeasurement.from_galcats([name_i+name_j]*2, psf_cat, q2_i, 'NZ_DUMMY', XX=gg, q2=q2_j, kernel_name2='NZ_DUMMY')
            psf_auto_specs.append(sp[0])
            varxi_arr+=list(varxis[0].copy())

    print 'varxi_arr',varxi_arr
    print 'cross_specs',cross_specs
    print 'psf_auto_specs',psf_auto_specs

    #make shape noise covariance
    shape_noise_cov = np.diag(np.array(varxi_arr))
    shape_noise_cov_info = twopoint.CovarianceMatrixInfo("COV_SHAPENOISE", [s.name for s in cross_specs]+[s.name for s in psf_auto_specs],
                                                         [len(s.value) for s in cross_specs]+[len(s.value) for s in psf_auto_specs], 
                                                         shape_noise_cov)
    #No covariance compuation, just save measurement
    t = twopoint.TwoPointFile(cross_specs + psf_auto_specs,  [twopoint.dummy_kernel("NZ_DUMMY")],
                                                      None, shape_noise_cov_info)
    t.to_fits(pj(args.outdir, 'corrs_covsn.fits'),clobber=True)
        

if __name__=="__main__":
    main()                 