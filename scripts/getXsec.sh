#!/bin/bash

myPWD=$PWD
logDir=$myPWD/logs/

###### the following variables are defined in a separate config file
# datasetFile=configs/datasets.dat
# jsonFile=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt
# jsonName=246908-260627-Prompt_25ns-golden_silver-v1
source configs/2015-v1.conf

datasets=(`cat $datasetFile | grep -v '#' | awk '{print $2}'`)
datasetNames=(`cat $datasetFile | grep -v '#' | awk '{print $1}'`)
IFS=$'\n'

#create temporary working area 
cd /tmp/$USER/
scram project CMSSW_7_5_0
cd CMSSW_7_5_0/src
eval `scramv1 runtime -sh`
wget -c https://raw.githubusercontent.com/syuvivida/generator/master/cross_section/runJob/ana.py  



if [ ! -d "tmp/" ];then mkdir tmp/; fi

for i in `seq 0 ${#datasets[@]}`
do
	dataset=${datasets[${i}]}
	datasetName=${datasetNames[${i}]}
	if [ -z "${datasetName}" ];then continue; fi # there is an empty line to be removed from the list

	echo $dataset
	l=`echo $line | awk '{print $0}'`
	case $dataset in
		*/MINIAODSIM)
			;;
		*)
			continue
			;;
	esac
	#echo $file
	mkdir -p $logDir/$dataset
	logfile=$logDir/$dataset/log.log
	if [ ! -e "$logfile" ];then
		file=`das_client --query="file dataset=$dataset" --limit=10 | grep -v "Showing"| sed 's|$|,|g'`
		file=`echo $file | sed -r 's|[,]+ |,|g;s|,,||;s|,$||'`
		#echo $file
		( echo $file &> $logfile; cmsRun ana.py inputFiles="$file" maxEvents=50000 &>> $logfile ) &
	fi

	
done


wait 

cd $myPWD

rm $logDir/l.list

for i in `seq 0 ${#datasets[@]}`
do
	dataset=${datasets[${i}]}
	datasetName=${datasetNames[${i}]}
	if [ -z "${datasetName}" ];then continue; fi # there is an empty line to be removed from the list


	case $dataset in
		*/MINIAODSIM)
			;;
		*)
			continue
			;;
	esac

	logfile=$logDir/$dataset/log.log

#Before Filtrer: total cross section = 1.015e+01 +- 2.284e-02 pb
#Filter efficiency (taking into account weights)= (50088) / (50088) = 1.000e+00 +- 0.000e+00
#Filter efficiency (event-level)= (50088) / (50088) = 1.000e+00 +- 0.000e+00
#After filter: final cross section = 1.015e+01 +- 2.284e-02 pb
	beforeFilter=`grep "Before matching" $logfile | cut -d '=' -f 2 | awk '{print $1, $3}'`
	filter_weight=`grep "Filter efficiency (taking into account weights)" $logfile | cut -d '=' -f 2 | awk '{print $1, $3}'`
	filter_event=`grep "Filter efficiency (event-level)" $logfile | cut -d '=' -f 2 | awk '{print $1, $3}'`
	afterFilter=`grep "After filter" $logfile | cut -d '=' -f 2 | awk '{print $1, "\t", $3}'`

	case $datasetName in
		*powheg*)
			afterFilter=$beforeFilter
			# echo $beforeFilter
			# echo $afterFilter
			;;
		*)
			;;
	esac
	echo $datasetName $afterFilter >> $logDir/l.list

	# this updated always
	#sed -i -r  "/$datasetName/{s|([[:alnum:]_]+\t[[:alnum:]_/-]+\t)[[:digit:].e+-]+\t[[:digit:].e+-]+\t|\1$afterFilter\t|}" $datasetFile

	# this updates only if not set
	sed -i -r  "/$datasetName/{s|([[:alnum:]_]+\t[[:alnum:]_/-]+\t)[-]\t[-]\t|\1$afterFilter\t|}" $datasetFile 
	./scripts/reformatDatasets.sh
done
