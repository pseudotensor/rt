#!/bin/sh
#BSUB -q kipac-ibq
#BSUB -u roman@astro.umd.edu
#BSUB -J qu1[34-34]
#BSUB -o transfer_%I.out
#BSUB -e transfer_%I.err
#BSUB -n 8
#BSUB -R "span[ptile=8]"

export OMP_NUM_THREADS=8
export LD_LIBRARY_PATH=/afs/slac.stanford.edu/package/intel_tools/2013u0/composer_xe_2013.0.079/compiler/lib/intel64
export PATH=/afs/slac.stanford.edu/package/intel_tools/2013u0/composer_xe_2013.0.079/compiler/lib/intel64

./ASTRORAY_main 1 3000 3049 1
# 1 - thickdisk7
# 3000 - minimum fieldline index
# 3049 - maximum fieldline index
# 1 - type of computation, quick computation
# LSB_JOBINDEX - initial conditions
# http://arxiv.org/abs/1007.4832
# bsub < qu1.bsub
