#!/bin/bash

hemin=0.0025
hemax=0.007
hestep=0.00025

himin=0.010
himax=0.025
histep=0.001



for he in $(seq $hemin $hestep $hemax); do
	for hi in $(seq $himin $histep $himax); do
	    #max=`echo $m + $pstep | bc`
	    program="python RUM_Detektor_HeHi_2ndversion_cluster.py -he $he -hi $hi"
	    echo $program        
	    
	    #bsub -q Batch24 $program
	done
done
