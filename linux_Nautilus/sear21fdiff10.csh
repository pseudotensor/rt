#PBS -A  TG-PHY120005
#PBS -l ncpus=16
#PBS -l mem=64000MB
#PBS -l walltime=23:59:00
cd /nics/d/home/shcher/rt

./transfer_quad 1 21 10 5
# 1 - thickdisk7
# 21 - number of shots to average over
# 10 - fdiff
# 5 - type of computation, search for minimum \chi^2
# PBS_ARRAYID - initial conditions
# http://arxiv.org/abs/1007.4832
# qsub -t 1-7 sear21.csh
