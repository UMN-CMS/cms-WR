#!/bin/bash

#all miniaod files (not skimmed) are listed in txt files in one directory
#no need to read skim file names from eos
#use the dataset tag defined in configs/ to label the dataset for later
#cross sxn times lumi weighting

#import datasetFile
source configs/2015-v1.conf

##for all datasets in configs/missing_datasets.dat or datasets.dat, depending on 2015-v1.conf
##first get all the dataset short names, and save them to a list named mcIdentifier
mcIdentifier=(`cat $datasetFile | grep -v '#' | awk '{print $1}'`)

#each list of files is stored in one dir in a txt file whose name starts with an element
#in mcIdentifier, and ends in Files.txt
inputFileListDir='/afs/cern.ch/work/s/skalafut/public/WR_starting2015/privateWRMiniAod/'
fileListEnd='Files.txt'

#now loop over all elements in mcIdentifier
for j in ${!mcIdentifier[*]}
do
	#number used to distinguish different jobs and input files tied to the same dataset
	startingCount=1

	##now get all the skim files from the dataset linked to mcIdentifier
	#inputFiles=`cat $inputFileListDir${mcIdentifier[$j]}$fileListEnd`

	##replace NNN by a number and INPTFILE with a file path name from inputFiles in temp_runAnalysis_cfg.py
	#eval "cd test"
	#eval "sed 's@NNN@$startingCount@g' temp_runAnalysis_cfg.py > tempOne.py"
	#eval "sed 's@TAGNAME@${mcIdentifier[$j]}@g' tempOne.py > tempTwo.py"
	#eval "sed 's@INPTFILE@$inputFiles@g' tempTwo.py > runAnalysis_cfg_${mcIdentifier[$j]}_${startingCount}.py"
	#rm tempOne.py tempTwo.py
	#eval "cd .."

	#replace NNN by a number, INPTFILELIST by an abs path name to a file in runAnalysisPrivateWR.sh
	eval "sed 's@NNN@$startingCount@g' runAnalysisPrivateWR.sh > tempOne.sh"
	eval "sed 's@TAGNAME@${mcIdentifier[$j]}@g' tempOne.sh > runAnalysisPrivateWR_${mcIdentifier[$j]}_${startingCount}.sh"
	eval "chmod u+x runAnalysisPrivateWR_${mcIdentifier[$j]}_${startingCount}.sh"
	rm tempOne.sh 

	#submit the job to 1nh queue and request at least 2 GB of disk
	#echo "bsub -R 'pool>2000' -q 1nh -J runAnalysisPrivateWR_${mcIdentifier[$j]}_Part_${startingCount} < runAnalysisPrivateWR_${mcIdentifier[$j]}_${startingCount}.sh"
	eval "bsub -R 'pool>4000' -q 8nh -J runAnalysisPrivateWR_${mcIdentifier[$j]}_Part_${startingCount} < runAnalysisPrivateWR_${mcIdentifier[$j]}_${startingCount}.sh"

	#let startingCount=startingCount+1

done

