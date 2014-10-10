#!/bin/bash
#PBS -A  TG-PHY120005
#PBS -l ncpus=88
#PBS -l mem=352000MB
#PBS -l walltime=23:59:00
cd /nics/d/home/shcher/rt

./ASTRORAY_main 1 1534 2033 1 &
sleep 5
./ASTRORAY_main 1 2034 2533 1 &
sleep 5
./ASTRORAY_main 1 2534 3033 1 &
sleep 5
./ASTRORAY_main 1 3034 3533 1 &
sleep 5
./ASTRORAY_main 1 3534 4033 1 &
sleep 5
./ASTRORAY_main 1 4034 4533 1 &
sleep 5
./ASTRORAY_main 1 4534 5033 1 &
sleep 5
./ASTRORAY_main 1 5034 5533 1 &
sleep 5
./ASTRORAY_main 1 5534 6033 1 &
sleep 5
./ASTRORAY_main 1 6034 6533 1 &
sleep 5
./ASTRORAY_main 1 6534 7034 1 &
wait

# 1 - thickdisk7
# 21 - number of shots to average over
# 2.0 - slope of B field
# 1 - type of computation, quick computation
# PBS_ARRAYID - change of magnetic field slope
# http://arxiv.org/abs/1007.4832
# qsub -t 28-28 qually.csh
