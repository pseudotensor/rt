#PBS -A  TG-PHY120005
#PBS -l ncpus=16
#PBS -l mem=8GB
#PBS -l walltime=00:40:00
cd /nics/d/home/shcher/rt

./ASTRORAY_main 1 2 2 5
# 1 - thickdisk7
# 2 - number of shots to average over
# 2 - fdiff
# 5 - type of computation, search for minimum \chi^2
# PBS_ARRAYID - initial conditions
# http://arxiv.org/abs/1007.4832
# qsub -t 1-1 sear2.csh
