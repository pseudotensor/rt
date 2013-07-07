#!/bin/bash
#PBS -A  TG-PHY120005
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=23:59:00
cd /nics/d/home/shcher/rt

./transfer_quad 1 3000 3332 4

# 1 - thickdisk7
# 21 - number of shots to average over
# 2.0 - slope of B field
# 4 - type of computation, imaging
# PBS_ARRAYID - change of magnetic field slope
# http://arxiv.org/abs/1007.4832
# qsub -t 30-32 imag31.csh
