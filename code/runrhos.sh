#!/bin/bash

echo "RUNNING"
echo "slurm procid = " $SLURM_PROCID
echo "slurm ntasks = " $SLURM_NTASKS


STARS=/global/cscratch1/sd/seccolf/y3a1-v29
OUT=/global/cscratch1/sd/alsina/alpha-beta-gamma/out

cmd="python reserved_stars_corr.py --piff_cat=$STARS --exps_file=ally3.grizY --outpath=$OUT --bands=riz"
date
echo $cmd
$cmd

