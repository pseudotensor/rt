#!/bin/bash
#PBS -A  TG-PHY120005
#PBS -l ncpus=80
#PBS -l mem=320000MB
#PBS -l walltime=23:59:00
cd /nics/d/home/shcher/rt

startdump=1435
enddump=6965
skipdump=70
dumpi=$startdump
while [ $dumpi -le $enddump ]
do
       dumpf=$(($dumpi+$skipdump-1))
	./transfer_quad 1 $dumpi $dumpf 1 &
       dumpi=$(($dumpf+1))
	sleep 10
done
wait

# 1 - thickdisk7
# 21 - number of shots to average over
# 2.0 - slope of B field
# 1 - type of computation, quick computation
# PBS_ARRAYID - change of magnetic field slope
# http://arxiv.org/abs/1007.4832
# qsub -t 28-28 quallz.csh
