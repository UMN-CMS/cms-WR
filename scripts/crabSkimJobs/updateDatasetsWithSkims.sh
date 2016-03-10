#!/bin/bash
if [ "$1" == "--force" ];then 
	FORCE=y
fi

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
	crabDir=crab/skim/crab_skim_${datasetName}${skimProductionTAG}/crab_$datasetName

	ls -d $crabDir &> /dev/null || 	crabDir=crab/skim/crab_skim_${datasetName}/crab_$datasetName

	#check if the skim of the sample have been submitted
	ls -d $crabDir > /dev/null || continue

	if [ ! "$(ls -A $crabDir/results)" ]; then 
		crab status -d $crabDir

		if [ "`grep -c COMPLETED $crabDir/crab.log`" != "0" -o "${FORCE}" == "y" ];then
			crab report -d $crabDir
		fi
	fi

	if [ "`grep -c COMPLETED $crabDir/crab.log`" != "0" -o "${FORCE}" == "y" ];then
		if [ ! -e "$crabDir/results/job_out.1.0.txt" ];then
			crab getlog --short --jobids=1-1000 --dir=$crabDir
		fi
		readEvents=`grep "events have been read" $crabDir/crab.log | cut -d ':' -f 4 | awk '{print $1}'`
		writtenEvents=`grep "events have been written" $crabDir/crab.log | cut -d ':' -f 4 | awk '{print $1}'`
		outputDataset=`grep "Output dataset:" $crabDir/crab.log | tail -1| cut -d ':' -f 2 | sed 's|[[:space:]]*||'`
		echo $datasetName $readEvents $writtenEvents $outputDataset
		sed -i -r  "/$datasetName/{s|([[:alnum:]_]+\t[[:alnum:]_/-]+\t[[:digit:].e+-]+\t[[:digit:].e+-]+\t)[[:digit:]-]+\t[[:alnum:]_/.+-]+\t[[:digit:].e+-]+|\1$readEvents\t$outputDataset\t$writtenEvents|}" $datasetFile 
		./scripts/reformatDatasets.sh
#		break
	fi
done
exit 0

for crabDir in crab_*/*
do
	
done
