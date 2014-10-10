#PBS -l walltime=04:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=4524mb
cd /a/f20-fs1/export/home/deepthought/roman/src

./ASTRORAY_main 1 21 2.0 2
# 1 - spin 0.
# 21 - number of shots to average over
# 2.0 - slope of B field
# 2 - type of computation, testing convergence
# PBS_ARRAYID - test number
# http://arxiv.org/abs/1007.4832
# qsub -A astroctc-hi -q narrow-extended -t 1-14 conv1.csh
