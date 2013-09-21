#!/bin/sh
#PBS -S /bin/sh
#PBS -A christoq_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -N SELF-STACK
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=40:00
#PBS -t 0-19
#PBS -V
#PBS -j oe
#

cd /nfs/christoq_ls/nkern/Documents/MDB_milliMil_halodata/Caustic

python 2D_gal_self_stack.py $PBS_ARRAYID 5 5 1 2






