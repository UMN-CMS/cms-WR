#!/bin/bash

#import productionTAG, skimProductionTAG and datasetFile
source configs/2015-v1.conf

##for all datasets in configs/missing_datasets.dat
##first get all the dataset short names, and save them to a vector
mcIdentifier=(`cat $datasetFile | grep -v '#' | awk '{print $1}'`)
fileTag='Files.txt'

#now loop over all elements in mcIdentifier
for j in ${!mcIdentifier[*]}
do
	#number used to distinguish different jobs processing the same dataset
	startingCount=1

	#now get all the files from the dataset linked to mcIdentifier
	fileList=(`cat ${mcIdentifier[$j]}$fileTag`)


	for i in ${!fileList[*]}
	do
		#replace NNN by a number and INPTFILE with a file path name from fileList
		eval "cd test"
		eval "sed 's@NNN@$startingCount@g' temp_skims_cfg.py > tempOne.py"
		eval "sed 's@TAGNAME@${mcIdentifier[$j]}@g' tempOne.py > tempTwo.py"
		eval "sed 's@INPTFILE@${fileList[$i]}@g' tempTwo.py > skims_cfg_${mcIdentifier[$j]}_${startingCount}.py"
		rm tempOne.py tempTwo.py
		eval "cd .."

		#update runSkimWJets.sh
		eval "sed 's@NNN@$startingCount@g' runSkimWJets.sh > tempOne.sh"
		eval "sed 's@TAGNAME@${mcIdentifier[$j]}@g' tempOne.sh > runSkimWJets_${mcIdentifier[$j]}_${startingCount}.sh"
		eval "chmod u+x runSkimWJets_${mcIdentifier[$j]}_${startingCount}.sh"
		rm tempOne.sh

		#submit the job to 1nh queue and request at least 2 GB of disk
		#echo "bsub -R 'pool>3000' -q 1nh -J runSkimWJets_${mcIdentifier[$j]}_Part_${startingCount} < runSkimWJets_${mcIdentifier[$j]}_${startingCount}.sh"
		eval "bsub -R 'pool>3000' -q 1nh -J runSkimWJets_${mcIdentifier[$j]}_Part_${startingCount} < runSkimWJets_${mcIdentifier[$j]}_${startingCount}.sh"

		let startingCount=startingCount+1
	done

done

