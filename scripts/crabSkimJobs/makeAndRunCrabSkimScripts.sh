#!/bin/bash
### this script should be done in python!!!! to be well integrated with the crab submission
#source /cvmfs/cms.cern.ch/crab3/crab.sh


###### the following variables are defined in a separate config file
# datasetFile=configs/datasets.dat
# jsonFile=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt
# jsonName=246908-260627-Prompt_25ns-golden_silver-v1
source configs/2015-v1.conf


crabFile=tmp/crab.py

datasets=(`cat $datasetFile | grep -v '#' | awk '{print $2}'`)
datasetNames=(`cat $datasetFile | grep -v '#' | awk '{print $1}'`)
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
config.General.workArea = 'crab/skim/crab_skim_'+"$datasetName${skimProductionTAG}"
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test/skims_cfg.py'
EOF

LUMI=2000
#### if the dataset is DATA or DY save the TagAndProbe triggers
case $datasetName in 
	DoubleEG*)
		params="$params, 'saveTnP=1', 'GT=74X_dataRun2_Prompt_v4'"
		LUMI=1000
		;;
	SingleMu*)
		params="$params, 'saveTnP=1', 'GT=74X_dataRun2_Prompt_v4'"
		LUMI=1000
		;;
	MuEG*)
		params="$params, 'saveTnP=0', 'GT=74X_dataRun2_Prompt_v4'"
		LUMI=1000
		;;
	DY*)
		params="$params, 'saveTnP=1', 'GT=74X_mcRun2_asymptotic_v4'"
		jsonFile=""
		;;
	
	*)
		params="$params, 'saveTnP=0', 'GT=74X_mcRun2_asymptotic_v4'"
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
config.Data.unitsPerJob = ${LUMI}

#True allows the jobs to run anywhere, regardless of where the input data is located
config.Data.ignoreLocality = False

#totalUnits only needs to be specified for GEN-SIM jobs
#config.Data.totalUnits = 200000
config.Data.outLFNDirBase = '/store/user/%s/WRskims/' % (getUsernameFromSiteDB())
config.Data.publication = True 
#config.Data.publishDataName = 'realData_FNLST_13TeV_CHNL_UNIQUE'
config.Data.outputDatasetTag =  config.General.requestName + "${skimProductionTAG}"
config.Data.lumiMask = "$jsonFile"


#a list of the only sites at which these jobs can run
#config.Site.whitelist = ["T2_US*"]
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T2_CH_CERN'

EOF


echo "crab submit -c $crabFile"

done
#skim_realData_${channel[$q]}_${datasets[$q]}_${suffix[$q]}.py"



