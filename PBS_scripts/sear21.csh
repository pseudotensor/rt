#PBS -A  TG-PHY120005
#PBS -l ncpus=16
#PBS -l mem=64GB
#PBS -l walltime=12:00:00
cd /nics/d/home/shcher/rt

./ASTRORAY_main 1 21 0 5
# 1 - thickdisk7
# 21 - number of shots to average over
# 0 - fdiff
# 5 - type of computation, search for minimum \chi^2
# PBS_ARRAYID - initial conditions
# http://arxiv.org/abs/1007.4832
# qsub -t 1-7 sear21.csh
