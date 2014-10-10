#!/bin/bash
st=70;
for startdump in {1365..6965..$st}
do
enddump=$(($startdump+$st-1))
echo $startdump > z.txt
echo $enddump > z.txt
done
