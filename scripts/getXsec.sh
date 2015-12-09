#!/bin/bash

myPWD=$PWD
datasetFile=configs/datasets.dat
datasets=`cat $datasetFile |grep -v '^#'`

logDir=$myPWD/logs/


#create temporary working area 
cd /tmp/$USER/
scram project CMSSW_7_5_0
cd CMSSW_7_5_0/src
eval `scramv1 runtime -sh`
curl https://raw.githubusercontent.com/syuvivida/generator/master/cross_section/runJob/ana.py  -o ana.py



IFS=$'\n'
for line in $datasets
do
	dataset=`echo $line | awk '{print $1}'`
	echo $dataset
	#echo $file
	mkdir -p $logDir/$dataset
	logfile=$logDir/$dataset/log.log
	if [ ! -e "$logfile" ];then
		file=`das_client --query="file dataset=$dataset" --limit=3 | tail -1 | sed 's|$|,|g'`
		file=`echo $file | sed 's|, |,|g; s|,$||'`
		( echo $file &> $logfile; cmsRun ana.py inputFiles="$file" maxEvents=50000 &>> $logfile ) &
	fi
	echo -n $dataset >> $logDir/l.list
	grep 'After filter' $logfile | cut -d '=' -f 2 | awk '{print "\t", $1, $3}' >> $logDir/l.list

done


wait 
