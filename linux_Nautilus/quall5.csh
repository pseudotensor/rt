#!/bin/bash
#PBS -A  TG-PHY120005
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=23:59:00
cd /nics/d/home/shcher/rt

./transfer_quad 1 5000 5999 1

# 1 - thickdisk7
# 21 - number of shots to average over
# 2.0 - slope of B field
# 1 - type of computation, quick computation
# PBS_ARRAYID - change of magnetic field slope
# http://arxiv.org/abs/1007.4832
# qsub -t 30-32 quall5.csh
