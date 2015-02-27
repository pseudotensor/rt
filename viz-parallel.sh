#####################################################################
## RUN THIS SCRIPT TO MAKE USE OF MULTIPLE CORES TO VISUALIZE DATA ##
## USAGE: ./viz-idev.sh [ITER_START] [ITER_END]                    ##
## DEPENDENCIES:                                                   ##
##               - ASSUMES ROMAN's PYTHON STUFF SITS IN ~/py/      ##
#####################################################################
#module load python ## LOAD PYTHON LIBRARIES

## POOR-MAN's PARALELLIZATION
export NR_OF_PROC=24 ## 4: dell-xps13 |16: STAMPEDE |24: bh01
export ITER_START=$1
export ITER_END=$2
export CASE_START=$3
export CASE_END=$4
TIMESTEPS_PER_CORE=$((($ITER_END-$ITER_START)/$NR_OF_PROC))

# GAP = 0
# if GAP;do 
## THIS ONE-LINER IS USEFUL FOR VISUALIZATION ONLY FEW MISSING SNAPSHOTS
# MISSING_FILES=`for i in $(seq -w 0 2601); do [ ! -f rho_xy$i.png ] && echo rho_xy$i.png; done`
# fi


## LOOP OVER DESIRED TIME STEPS ##
for i in `seq $ITER_START 1 $ITER_END`; do
  for case in `seq $CASE_START 1 $CASE_END`;do 

    #filename=echo "rho_xy`seq -f "%04g" $i $i`.png"
    #filename=`echo "b2-over-rho-"\`seq -f "%04g" $i $i\`".png"`
    #if [ -f $filename ]; then 
	#echo "SKIPPING ALREADY EXISTING FILE AT ITER "$i
    #else
        # VISUALIZE, SEND TO BACKGROUND, TIME IT (NO ";" AT END OF LINE!)
	/usr/bin/time -f "SNAPSHOT $i CASE $case TOOK %e SECS" python rt-imaging.py shotimag93.75th140f230fn$i\_299.dat vary t & # vary t # vary r # vary theta # vary magn_cap # vary magn_floor
    #fi

    # WHEN RUNNING OUT OF PROCESSORS WAIT FOR BACKGROUND JOBS TO FINISH
    if [[ $(( ($i-$ITER_START)%$NR_OF_PROC )) -eq $(($NR_OF_PROC-1)) ]];then wait;fi

  done #
done ##
######

wait # ... UNTIL RETURN THE CMD LINE PROMPT
echo "|== DONE ==|"
