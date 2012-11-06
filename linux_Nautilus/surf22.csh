#PBS -l walltime=32:00:00
#PBS -l nodes=1:ppn=4
cd /a/f20-fs1/export/home/deepthought/roman/src

./transfer_quad 2 2 2.0 3
# 2 - spin 0.5
# 2 - type of surfing: accretion rate and heating constant
# 2.0 - slope of B field
# 3 - type of computation: surfing near best fit
# PBS_ARRAYID - number of snapshot: 0..10
# http://arxiv.org/abs/1007.4832
# qsub -A astroctc-hi -q narrow-long -t 0-10 surf22.csh
