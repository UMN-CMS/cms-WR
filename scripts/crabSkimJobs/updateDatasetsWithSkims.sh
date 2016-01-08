#!/bin/bash

###### the following variables are defined in a separate config file
# datasetFile=configs/datasets.dat
# jsonFile=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt
# jsonName=246908-260627-Prompt_25ns-golden_silver-v1
source configs/2015-v1.conf

datasets=(`cat $datasetFile | grep -v '#' | awk '{print $2}'`)
datasetNames=(`cat $datasetFile | grep -v '#' | awk '(NF==2 || NF>=4){print $1}'`)
IFS=$'\n'

if [ ! -d "tmp/" ];then mkdir tmp/; fi

for i in `seq 0 ${#datasets[@]}`
do
	dataset=${datasets[${i}]}
	datasetName=${datasetNames[${i}]}
	if [ -z "${datasetName}" ];then continue; fi # there is an empty line to be removed from the list
	crabDir=crab_skim_$datasetName/crab_$datasetName

	#check if the skim of the sample have been submitted
	ls -d $crabDir > /dev/null || continue

	if [ ! "$(ls -A $crabDir/results)" ]; then 
		crab status -d $crabDir
		if [ "`grep -c COMPLETED $crabDir/crab.log`" != "0" ];then
			crab report -d $crabDir
		fi
	fi

	if [ "`grep -c COMPLETED $crabDir/crab.log`" != "0" ];then
		readEvents=`grep "events have been read" $crabDir/crab.log | cut -d ':' -f 4 | awk '{print $1}'`
		echo $datasetName $readEvents
		sed -i -r  "/$datasetName/{s|([[:alnum:]_]+\t[[:alnum:]_/-]+\t[[:digit:].-]+\t[[:digit:].-]+\t)[[:digit:].-]+|\1$readEvents\t|}" $datasetFile 
		./scripts/reformatDatasets.sh
#		break
	fi
done
exit 0

for crabDir in crab_*/*
do
	
done