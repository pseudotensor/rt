#!/bin/bash
#PBS -A  TG-PHY120005
#PBS -l ncpus=16
#PBS -l mem=16GB
#PBS -l walltime=00:59:00
cd /nics/d/home/shcher/rt

./transfer_quad 1 1965 1999 1

# 1 - thickdisk7
# 21 - number of shots to average over
# 2.0 - slope of B field
# 1 - type of computation, quick computation
# PBS_ARRAYID - change of magnetic field slope
# http://arxiv.org/abs/1007.4832
# qsub -t 28-28 quall1.csh
