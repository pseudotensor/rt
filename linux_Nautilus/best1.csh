#PBS -l walltime=5:00:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=4524mb
cd /a/f20-fs1/export/home/deepthought/roman/src

./transfer_quad 1 21 2.0 5
# 1 - spin 0.
# 21 - number of shots to average over
# 2.0 - slope of B field
# 5 - type of computation, search for the best fit
# PBS_ARRAYID - initial guess of the best parameters
# http://arxiv.org/abs/1007.4832
# qsub -A astroctc-hi -q narrow-extended -t 1-7 best1.csh
