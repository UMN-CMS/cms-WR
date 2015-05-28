#!/bin/bash

nuMasses=(50 100 200 300)
wrMasses=(2600 2600 2600 2600)
avgPileup=40
inputDatasets=('/testDatasetOne/USER' '/testDatasetTwo/USER')

for h in ${!inputDatasets[*]}
do
	echo 'crab submit -c WR_MWr${wrMasses[$h]}_to_ENu_MNu${nuMasses[$h]}_13TeV_DIGIL1DIGI2RAW_RAW2DIGIL1RecoRECO_${avgPileup}PU_crab.py'
	#eval 'crab submit -c WR_MWr${wrMasses[$h]}_to_ENu_MNu${nuMasses[$h]}_13TeV_DIGIL1DIGI2RAW_RAW2DIGIL1RecoRECO_${avgPileup}PU_crab.py'

done

