#!/bin/bash

#use integer values for masses and pileup
nuMasses=(50 100 200 300)
wrMasses=(2600 2600 2600 2600)
avgPileup=40
inputDatasets=('/testDatasetOne/USER' '/testDatasetTwo/USER')

#echo "number of input datasets: ${#inputDatasets[*]}"
#echo "length of first dataset name: ${#inputDatasets[0]}"
for h in ${!inputDatasets[*]}
do
	#make crab config file
	eval "sed 's/MWrUndef/MWr${wrMasses[$h]}/g' WR_MWrUndef_to_ENu_MNuUndef_13TeV_DIGIL1DIGI2RAW_RAW2DIGIL1RecoRECO_NNUUMMPU_crab.py > WR_tempOne_crab.py"
	eval "sed 's/MNuUndef/MNu${nuMasses[$h]}/g'  WR_tempOne_crab.py > WR_tempTwo_crab.py"
	eval "sed 's/NNUUMM/$avgPileup/g' WR_tempTwo_crab.py > WR_tempThree_crab.py"
	
	#the delimiter used in the next sed command cannot be a forward slash
	eval "sed 's@datasetFromDBS@${inputDatasets[$h]}@g' WR_tempThree_crab.py > WR_MWr${wrMasses[$h]}_to_ENu_MNu${nuMasses[$h]}_13TeV_DIGIL1DIGI2RAW_RAW2DIGIL1RecoRECO_${avgPileup}PU_crab.py"
	rm WR_tempOne_crab.py WR_tempTwo_crab.py WR_tempThree_crab.py

	#make python config file which governs the DIGI,L1,RAW,L1Reco,and RECO simulation 
	eval "sed 's/MWrUndef/MWr${wrMasses[$h]}/g' WR_MWrUndef_ToENu_MNuUndef_13TeV_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_NNUUMMPU.py > WR_tempOne.py"
	eval "sed 's/MNuUndef/MNu${nuMasses[$h]}/g'  WR_tempOne.py > WR_tempTwo.py"
	eval "sed 's/NNUUMM/$avgPileup/g'  WR_tempTwo.py > WR_MWr${wrMasses[$h]}_ToENu_MNu${nuMasses[$h]}_13TeV_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_${avgPileup}PU.py"
	rm WR_tempOne.py WR_tempTwo.py

done

