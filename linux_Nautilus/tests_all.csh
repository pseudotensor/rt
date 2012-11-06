#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=4524mb
cd /a/f20-fs1/export/home/deepthought/roman/src

./transfer_quad 2 21 2.0 1
# 2 - spin 0.5
# 21 - number of shots to average over
# 2.0 - slope of B field
# 1 - type of computation, quick computation w/ effects switched on/off
# http://arxiv.org/abs/1007.4832
# qsub -A astroctc-hi -q narrow-med -t 7-11 tests_all.csh
