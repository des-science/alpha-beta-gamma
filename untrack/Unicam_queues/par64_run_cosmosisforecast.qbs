#!/bin/bash
#PBS -M a.navarro.alsina@gmail.com
#PBS -m abe
#PBS -N par64a
#PBS -e par64a.err
#PBS -o par64a.out
#PBS -q par64
#PBS -l nodes=4:ppn=16

source /home/sw/masternode/intel/2015/install/composerxe/bin/compilervars.sh intel64
source /home/sw/masternode/intel/2015/install/mpi/impi/5.1.2.150/bin64/mpivars.sh

export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/opt/pbs/bin/pbs_tmrsh
export I_MPI_DEVICE=rdssm

export PATH=/home/dfa/sobreira/alsina/sw/pyhton/2714/install/bin:$PATH
export PYTHONPATH=$PYTHONPATH:/home/dfa/sobreira/alsina/sw/galsim/install/lib/python2.7/site-packages
export PATH=/home/dfa/sobreira/alsina/sw/cfitsio/install/bin:$PATH
export LD_LIBRARY_PATH=/home/dfa/sobreira/alsina/sw/cfitsio/install/lib:/home/dfa/sobreira/alsina/sw/ccfits/25/install/lib:/home/dfa/sobreira/alsina/sw/tmv/install/lib:/home/dfa/sobreira/alsina/sw/boost/166/install/lib:$LD_LIBRARY_PATH

cd /home/dfa/sobreira/alsina/sw/cosmosis/cosmosis/wl-pipeline-andres-debug
source /home/dfa/sobreira/alsina/sw/cosmosis/cosmosis/setup-my-cosmosis

INSTALL=/home/dfa/sobreira/alsina/sw
START_PATH=/home/dfa/sobreira/alsina/alpha-beta-gamma/cosmosis_pipe

mpirun --mca btl vader,tcp,self -n 16 $INSTALL/cosmosis/cosmosis/bin/cosmosis --mpi $START_PATH/run_d_l.ini
