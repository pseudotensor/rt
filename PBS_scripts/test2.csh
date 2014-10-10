#!/bin/bash
startdump=1435
enddump=6965
skipdump=70
dumpi=$startdump
while [ $dumpi -le $enddump ]
do
       dumpf=$(($dumpi+$skipdump-1))
	echo $dumpi
	echo $dumpf
       dumpi=$(($dumpf+1))
done
