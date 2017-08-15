#!/bin/bash

hemin=0.001
hemax=0.010
hestep=0.0005

himin=0.015
himax=0.025
histep=0.001



for he in $(seq $hemin $hestep $hemax); do
	for hi in $(seq $himin $histep $himax); do
	    #max=`echo $m + $pstep | bc`
	    program="python RUM_Detektor_2Nii.py -he $he -hi $hi"
	    echo $program        
	    
	    bsub -q Batch24 $program
	done
done
