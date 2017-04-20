#!/bin/bash

#import productionTAG, skimProductionTAG and datasetFile
source configs/2015-v1.conf

##for all datasets in configs/missing_datasets.dat or datasets.dat, depending on 2015-v1.conf
##first get all the dataset short names, and save them to a list named mcIdentifier
mcIdentifier=(`cat $datasetFile | grep -v '#' | awk '{print $1}'`)
eosReadingTag='root://eoscms.cern.ch/'
extraTag='_NoLeptId'

#path to skim files
eosSkimPath='/eos/cms/store/user/shervin/skims/'

#now loop over all elements in mcIdentifier
for j in ${!mcIdentifier[*]}
do
	#number used to distinguish different jobs processing the same dataset
	startingCount=1

	#now get all the skim files from the dataset linked to mcIdentifier
	inputFiles=(`eos ls $eosSkimPath${mcIdentifier[$j]}$skimProductionTAG/`)

	for i in ${!inputFiles[*]}
	do
		#echo $eosSkimPath${mcIdentifier[$j]}$skimProductionTAG/${inputFiles[$i]}
		#replace NNN by a number and INPTFILE with a file path name from inputFiles in temp_runAnalysis_cfg.py
		eval "cd test"
		eval "sed 's@NNN@$startingCount@g' temp_runAnalysis_cfg.py > tempOne.py"
		eval "sed 's@TAGNAME@${mcIdentifier[$j]}$extraTag@g' tempOne.py > tempTwo.py"
		eval "sed 's@INPTFILE@$eosReadingTag$eosSkimPath${mcIdentifier[$j]}$skimProductionTAG/${inputFiles[$i]}@g' tempTwo.py > runAnalysis_cfg_${mcIdentifier[$j]}_${startingCount}.py"
		rm tempOne.py tempTwo.py
		eval "cd .."

		#replace NNN by a number in runAnalysisQcdFromData.sh
		eval "sed 's@NNN@$startingCount@g' runAnalysisQcdFromData.sh > tempOne.sh"
		eval "sed 's@TAGNAME@${mcIdentifier[$j]}@g' tempOne.sh > runAnalysisQcdFromData_${mcIdentifier[$j]}_${startingCount}.sh"
		eval "chmod u+x runAnalysisQcdFromData_${mcIdentifier[$j]}_${startingCount}.sh"
		rm tempOne.sh

		#submit the job to 1nh queue and request at least 2 GB of disk
		#echo "bsub -R 'pool>2000' -q 1nh -J runAnalysisQcdFromData_${mcIdentifier[$j]}_Part_${startingCount} < runAnalysisQcdFromData_${mcIdentifier[$j]}_${startingCount}.sh"
		eval "bsub -R 'pool>2000' -q 1nh -J runAnalysisQcdFromData_${mcIdentifier[$j]}_Part_${startingCount} < runAnalysisQcdFromData_${mcIdentifier[$j]}_${startingCount}.sh"

		let startingCount=startingCount+1
	done

done

