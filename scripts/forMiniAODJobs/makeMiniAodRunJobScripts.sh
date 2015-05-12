#!/bin/bash

for q in {2..26} 
do
	#replace ABCDE and NUM in .csh file and reco_and_miniAOD.py file
	#replace NUM in job submission file

	eval "sed 's/NUM/$q/g' runMiniAODMaker.csh > runMiniAODMaker_temp.csh"
	eval "sed 's/ABCDE/WR_signal/g' runMiniAODMaker_temp.csh > runMiniAODMaker_$q.csh"
	
	eval "sed 's/NUM/$q/g' miniAODMaker > miniAODMaker_$q"
	
	eval "sed 's/NUM/$q/g' reco_and_miniAOD.py > reco_and_miniAOD_temp.py"
	eval "sed 's/ABCDE/WR_signal/g' reco_and_miniAOD_temp.py > reco_and_miniAOD_$q.py"

	rm runMiniAODMaker_temp.csh reco_and_miniAOD_temp.py 

done


