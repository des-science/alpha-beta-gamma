#!/bin/bash

echo "RUNNING"
echo "slurm procid = " $SLURM_PROCID
echo "slurm ntasks = " $SLURM_NTASKS


STARS=/global/cscratch1/sd/seccolf/y3a1-v29
METACAL=/global/cscratch1/sd/troxel/cats_des_y3/Y3_mastercat_v2_6_20_18.h5
OUT=/global/cscratch1/sd/alsina/alpha-beta-gamma/out

cmd="python galaxies_rstars_corr.py --metacal_cat=$METACAL --piff_cat=$STARS --exps_file=ally3.grizY --outpath=$OUT --bands=riz"
date
echo $cmd
$cmd



