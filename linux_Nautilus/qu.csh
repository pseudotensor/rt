#PBS -l nodes=1:ppn=4
cd /a/f20-fs1/export/home/deepthought/roman/src

./transfer_quad 2 2 2.0 1
# 2 - spin 0.5
# 2 - number of shots to average over
# 2.0 - slope of B field
# 1 - type of computation, quick computation
# PBS_ARRAYID - change of magnetic field slope
# http://arxiv.org/abs/1007.4832
# qsub -A astroctc-hi -q debug -t 0 qu.csh
