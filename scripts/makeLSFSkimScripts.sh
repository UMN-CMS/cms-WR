#!/bin/bash
### this script should be done in python!!!! to be well integrated with the crab submission
#source /cvmfs/cms.cern.ch/crab3/crab.sh
CREATE=y
SUBMIT=y
SCHEDULER=caf
#SCHEDULER=remoteGlidein

###### the following variables are defined in a separate config file
# datasetFile=configs/datasets.dat
# jsonFile=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt
# jsonName=246908-260627-Prompt_25ns-golden_silver-v1
source configs/2015-v1.conf


setStoragePath(){
    #$1 = storage_element
    #$2 = scheduler
    #echo "[DEBUG] Setting storage path for $1 with scheduler $2"
    case $1 in
        caf* | T2_CH_CERN)
            case $2 in
                caf|lsf)
                    STORAGE_ELEMENT=caf.cern.ch
                    STORAGE_PATH=root://eoscms//eos/cms/store
                    ;;
            #glite | glidein)
                remoteGlidein|condor)
                    STORAGE_ELEMENT=srm-eoscms.cern.ch
                    STORAGE_PATH=/srm/v2/server?SFN=/eos/cms/store
                #STORAGE_ELEMENT=caf.cern.ch
                #STORAGE_PATH=root://eoscms//eos/cms/store
                    ;;
            *)
                    echo "[ERROR] Scheduler $2 for storage_element $1 not implemented" >> /dev/stderr
                    exit 1
                    ;;
            esac
            ;;
        T2_IT_Rome)
            STORAGE_PATH=/srm/managerv2?SFN=/castor/cern.ch
            ;;
        clusterCERN)
            STORAGE_PATH=root://pccmsrm27.cern.ch:1094//cms/local
            ;;
    esac
}


#------------------------------ parsing

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -u -o h -l help,ui_working_dir:,createOnly,submitOnly,check,scheduler:,datasetName: -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    exit 1
fi


set -- $options
#echo $options

while [ $# -gt 0 ]
do
    case $1 in
		-h|--help) usage; exit 0;;
		--createOnly) echo "[OPTION] createOnly"; unset SUBMIT;;
		--submitOnly) echo "[OPTION] submitOnly"; unset CREATE;;
		--scheduler)  echo "[OPTION] scheduler = $2"; SCHEDULER=$2; shift;;
		--datasetName) echo "[OPTION] dataset = $2"; DATASETNAME=$2; shift;;
		--check) CHECK=y; EXTRAOPTION="--check"; unset CREATE; unset SUBMIT;;
		(--) shift; break;;
		(-*) usage; echo "$0: error - unrecognized option $1" 1>&2; usage >> /dev/stderr; exit 1;;
		(*) break;;
    esac
    shift
done


if [ ! -d "crab/analysis" ]; then mkdir crab/analysis -p; fi

crabFile=tmp/crab.py
crab2File=tmp/crab2.cfg
#DBS_URL=global

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


	if [ -z "${datasetName}" -o "${dataset}" == "-" ];then continue; fi
	if [ -n "${DATASETNAME}" -a "${datasetName}" != "${DATASETNAME}" ];then continue; fi

	UI_WORKING_DIR=crab/skims/crab_skims_${datasetName}${skimProductionTAG}

	setStoragePath caf $SCHEDULER

#	STORAGE_ELEMENT=caf.cern.ch
#	STORAGE_PATH=root://eoscms//eos/cms/store
	USER_REMOTE_DIR=/user/shervin/skims/${datasetName}${skimProductionTAG}/

	OUTFILES=${datasetName}.root

	if [ -n "${CREATE}" ];then
#		echo $dataset $datasetName
#		echo $i $datasetName
	#echo $dataset
		crabFile=tmp/skim_$datasetName.py
		crab2File=tmp/skim_$datasetName.cfg

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
pars=`echo $params | sed 's|,||g;' | tr -d \' `
echo $pars
#exit 0

		cat > $crab2File <<EOF
[CRAB]
jobtype = cmssw
scheduler = $SCHEDULER

[LSF]
queue = 1nd
[CAF]
queue = cmscaf1nd
resource = type==SLC6_64 


[CMSSW]
#allow_NonProductionCMSSW = 1

pset=test/skims_cfg.py
pycfg_params=$pars output=${OUTFILES}

#runselection=${RUNRANGE}
split_by_run=0

#output_file=${OUTFILES}
get_edm_output=1
check_user_remote_dir=1

datasetpath=${dataset}
EOF

		case $isMC in
			1)
				cat >> $crab2File <<EOF
total_number_of_events = -1
events_per_job=50000
EOF
				;;
			*)
				cat >> $crab2File <<EOF
total_number_of_lumis = -1
lumis_per_job=${LUMI}
EOF
				;;
		esac

		cat >> $crab2File <<EOF
#dbs_url = ${DBS_URL}

[USER]
ui_working_dir=$UI_WORKING_DIR
return_data = 0
copy_data = 1

storage_element=${STORAGE_ELEMENT}
user_remote_dir=$USER_REMOTE_DIR
storage_path=$STORAGE_PATH

thresholdLevel=50
eMail = shervin@cern.ch

[GRID]

rb = HC
rb = CERN
proxy_server = myproxy.cern.ch

EOF


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



#echo "done: $crab2File created"



crab -cfg ${crab2File} -create || {
	rm ${UI_WORKING_DIR} -Rf
	$0 --scheduler=remoteGlidein --datasetName=${datasetName}
}	
#exit 1
fi
if [ -n "${SUBMIT}" ];then
	crab -c ${UI_WORKING_DIR} -submit #|| exit 1
fi

if [ -n "${CHECK}" ];then
    if [ ! -e "${UI_WORKING_DIR}/res/finished" ];then
		echo "[STATUS] Unfinished ${UI_WORKING_DIR}"
		resubmitCrab.sh -u ${UI_WORKING_DIR}
    else
		for file in $OUTFILES
		do
			file=`basename $file .root`
			#mergeOutput.sh -u ${UI_WORKING_DIR} -g $file
		done
		
#		echo "${datasetName}" 
		if [ ! -e "${UI_WORKING_DIR}/res/lumiSummary.json" ];then 
			crab -c ${UI_WORKING_DIR} -report
        readEvents=`grep "Total Events read:" ${UI_WORKING_DIR}/log/crab.log | tail -1 | cut -d ':' -f 2 | awk '{print $1}'`
		writtenEvents="-"
#		writtenEvents=`grep "events have been written" $crabDir/crab.log | cut -d ':' -f 4 | awk '{print $1}'`
        outputDataset=`echo $STORAGE_PATH/${USER_REMOTE_DIR} | sed 's|/srm/v2/server?SFN=||'`
#		outputDataset=`grep "Output dataset:" $crabDir/crab.log | tail -1| cut -d ':' -f 2 | sed 's|[[:space:]]*||'`
#		echo $datasetName $readEvents $writtenEvents $outputDataset
		sed -i -r  "/$datasetName/{s|([[:alnum:]_]+\t[[:alnum:]_/-]+\t[[:digit:].e+-]+\t[[:digit:].e+-]+\t)[[:digit:]-]+\t[[:alnum:]_/.+-]+\t[[:digit:].e+-]+|\1$readEvents\t$outputDataset\t$writtenEvents|}" $datasetFile 
		./scripts/reformatDatasets.sh
		fi
    fi
#    echo "mergeOutput.sh -u ${UI_WORKING_DIR} -n ${DATASETNAME} -r ${RUNRANGE}"
fi


done




