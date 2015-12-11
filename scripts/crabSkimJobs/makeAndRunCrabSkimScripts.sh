#!/bin/bash
### this script should be done in python!!!! to be well integrated with the crab submission
#source /cvmfs/cms.cern.ch/crab3/crab.sh

datasetFile=configs/datasets.dat
jsonFile=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt
jsonName=246908-260627-Prompt_25ns-golden_silver-v1
crabFile=tmp/crab.py


datasets=(`cat $datasetFile | grep -v '#' | awk '{print $1}'`)
datasetNames=(`cat $datasetFile | grep -v '#' | awk '(NF==2 || NF>=4){print $2}'`)
IFS=$'\n'

if [ ! -d "tmp/" ];then mkdir tmp/; fi

#echo ${#datasets[*]}
for i in `seq 0 ${#datasets[@]}`
do
	params=""
	dataset=${datasets[${i}]}
	datasetName=${datasetNames[${i}]}
	if [ -z "${datasetName}" ];then continue; fi
	echo $i $datasetName
	#echo $dataset
	crabFile=tmp/skim_$datasetName.py

	cat > $crabFile  <<EOF
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = "$datasetName"
config.General.workArea = 'crab_skim_'+"$datasetName"
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test/skims_cfg.py'
EOF

#### if the dataset is DATA or DY save the TagAndProbe triggers
case $dataset in 
	/DoubleEG/*)
		params="$params, 'saveTnP=1', 'GT=74X_dataRun2_Prompt_v4'"
		;;
	/DoubleMu/*)
		params="$params, 'saveTnP=1', 'GT=74X_dataRun2_Prompt_v4'"
		;;
	/MuonEG/*)
		params="$params, 'saveTnP=1', 'GT=74X_dataRun2_Prompt_v4'"
		;;
	DY*)
		params="$params, 'saveTnP=1', 'GT=74X_mcRun2_asymptotic_v2'"
		jsonFile=""
		;;
	*)
		params="$params, 'saveTnP=0', 'GT=74X_mcRun2_asymptotic_v2'"
		jsonFile=""
		;;
esac
params=`echo $params | sed -r 's|^,||;s|[,]+|,|g'`

cat >> $crabFile <<EOF
config.JobType.pyCfgParams = [ $params ]
#config.JobType.maxMemoryMB = 2500 should not need this option for skims

config.Data.inputDataset = "$dataset"
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'LUMI'
config.Data.unitsPerJob = 30 

#True allows the jobs to run anywhere, regardless of where the input data is located
config.Data.ignoreLocality = True

#totalUnits only needs to be specified for GEN-SIM jobs
#config.Data.totalUnits = 200000
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True 
#config.Data.publishDataName = 'realData_FNLST_13TeV_CHNL_UNIQUE'
config.Data.outputDatasetTag =  config.General.requestName
config.Data.lumiMask = "$jsonFile"


#a list of the only sites at which these jobs can run
#config.Site.whitelist = ["T2_US*"]
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T2_CH_CERN'

EOF


echo "crab submit -c $crabFile"

done
#skim_realData_${channel[$q]}_${datasets[$q]}_${suffix[$q]}.py"



