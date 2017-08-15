#!/bin/bash

hemin=0.001
hemax=0.007
hestep=0.0005

himin=0.002
himax=0.008
histep=0.0005



for he in $(seq $hemin $hestep $hemax); do
	for hi in $(seq $himin $histep $himax); do
	    #max=`echo $m + $pstep | bc`
	    program="python RUM_Detektor_2Npp.py -he $he -hi $hi"
	    echo $program        
	    
	    bsub -q Batch24 $program
	done
done
