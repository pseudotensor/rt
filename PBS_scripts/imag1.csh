#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=4524mb
cd /a/f20-fs1/export/home/deepthought/roman/src

./ASTRORAY_main 1 0 2.0 4
# 1 - spin 0
# 0 - not used
# 2.0 - slope of B field
# 4 - type of computation: imaging
# PBS_ARRAYID - snapshot number
# http://arxiv.org/abs/1007.4832
# qsub -A astroctc-hi -q narrow-med -t 0-20 imag1.csh
