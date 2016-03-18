#!/bin/bash
### this script should be done in python!!!! to be well integrated with the crab submission
#source /cvmfs/cms.cern.ch/crab3/crab.sh
CREATE=y
SUBMIT=y
FILE_PER_JOB=2

###### the following variables are defined in a separate config file
# datasetFile=configs/datasets.dat
# jsonFile=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt
# jsonName=246908-260627-Prompt_25ns-golden_silver-v1
source configs/2015-v1.conf



#------------------------------ parsing

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -u -o h -l help,ui_working_dir:,createOnly,submitOnly,check,datasetName: -- "$@")
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
		--check) CHECK=y; EXTRAOPTION="--check"; unset CREATE; unset SUBMIT;;
		--datasetName) echo "[OPTION] dataset = $2"; DATASETNAME=$2; shift;;
		(--) shift; break;;
		(-*) usage; echo "$0: error - unrecognized option $1" 1>&2; usage >> /dev/stderr; exit 1;;
		(*) break;;
    esac
    shift
done


if [ ! -d "crab/analysis" ]; then mkdir crab/analysis -p; fi

crabFile=tmp/crab.py
crab2File=tmp/crab2.cfg
DBS_URL=phys03

datasets=(`cat $datasetFile | grep -v '#' | awk '{print $6}'`) # these are the skimmed datasets!
datasetNames=(`cat $datasetFile | grep -v '#' | awk '{print $1}'`)
IFS=$'\n'

if [ ! -d "tmp/" ];then mkdir tmp/; fi

#echo ${#datasets[*]}
for i in `seq 0 ${#datasets[@]}`
do
	params=""
	dataset=${datasets[${i}]}
	datasetName=${datasetNames[${i}]}
	if [ -n "${DATASETNAME}" -a "${datasetName}" != "${DATASETNAME}" ];then continue; fi



	if [ -z "${datasetName}" -o "${dataset}" == "-" ];then continue; fi
	if [ "`echo ${dataset} |grep -c ${skimProductionTAG}`" == "0" ];then 
		echo "skipping $datasetName because old skimProductionTAG" 
		continue 
	fi

	UI_WORKING_DIR=crab/analysis/crab_analysis_${datasetName}${productionTAG}
	
	STORAGE_ELEMENT=caf.cern.ch
	STORAGE_PATH=root://eoscms//eos/cms/store
	USER_REMOTE_DIR=/user/shervin/ntuples/${datasetName}${productionTAG}/unmerged
	SCHEDULER=caf

	OUTFILES=${datasetName}.root

	if [ -n "${CREATE}" -a ! -d "${UI_WORKING_DIR}" ];then
#		echo $dataset $datasetName
#		echo $i $datasetName
	#echo $dataset
		crabFile=tmp/analysis_$datasetName.py
		crab2File=tmp/analysis_$datasetName.cfg

		LUMI=2000
#### if the dataset is DATA or DY save the TagAndProbe triggers
		case $dataset in 
			/DoubleEG/*)
				params="$params, 'isMC=0'"
				isMC=0
				;;
			/SingleMuon/*)
				params="$params, 'isMC=0'"
				isMC=0
				;;
			/MuonEG/*)
				params="$params, 'isMC=0'"
				isMC=0
				;;
			/DY*)
				params="$params, 'isMC=1'"
				jsonFile=""
				isMC=1
				;;
			
			*)
				params="$params, 'isMC=1'"
				isMC=1
				jsonFile=""
				;;
		esac
		params=`echo $params | sed -r 's|^,||;s|[,]+|,|g'`



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

pset=test/runAnalysis_cfg.py
pycfg_params=isMC=${isMC} datasetTag=${datasetName}

#runselection=${RUNRANGE}
split_by_run=0

#output_file=${OUTFILES}
get_edm_output=0
check_user_remote_dir=1
EOF
		case ${dataset} in
			*/eos/*)
				echo ${datasetName} $dataset
				makefilelist.sh "${datasetName}"  ${dataset}  || exit 1

				FILELIST=filelist/${datasetName}.list
				if [ -n "$FILELIST" ]; then
					nFiles=`cat $FILELIST | wc -l`
					if [ -n "$NJOBS" ];then
						let FILE_PER_JOB=$nFiles/$NJOBS
						if [ "`echo \"$nFiles%$NJOBS\" | bc`" != "0" ];then
							let FILE_PER_JOB=$FILE_PER_JOB+1
						fi
					elif [ -n "$FILE_PER_JOB" ];then
						NJOBS=`perl -w -e "use POSIX; print ceil($nFiles/${FILE_PER_JOB}), qq{\n}"`
						if [ "`echo \"${nFiles}%${FILE_PER_JOB}\" | bc -l`" != "0" ];then
							let NJOBS=$NJOBS+1
						fi
					else
						NJOBS=$nFiles
						FILE_PER_JOB=1
					fi
				fi

				cat >> $crab2File <<EOF
datasetpath=None
total_number_of_events=${NJOBS}
number_of_jobs=${NJOBS}
EOF

				if [ "${isMC}" != "1" ];then
					cat >> ${crab2File} <<EOF
lumi_mask=${jsonFile}
EOF
				fi
					;;
			*)
				cat >> ${crab2File} <<EOF
dbs_url = ${DBS_URL}
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
				;;
		esac

		cat >> $crab2File <<EOF


[USER]
ui_working_dir=$UI_WORKING_DIR
return_data = 0
copy_data = 1

storage_element=caf.cern.ch
user_remote_dir=$USER_REMOTE_DIR
storage_path=$STORAGE_PATH

thresholdLevel=50
eMail = shervin@cern.ch

[GRID]

rb = HC
rb = CERN
proxy_server = myproxy.cern.ch

EOF


		cat > $crabFile  <<EOF
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = "$datasetName"
config.General.workArea = 'crab/analysis/crab_analysis_'+"$datasetName"
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test/runAnalysis_cfg.py'
EOF


		cat >> $crabFile <<EOF
config.JobType.pyCfgParams = [ $params, "datasetTag=$datasetName" ]
#config.JobType.maxMemoryMB = 2500 should not need this option for skims

config.Data.inputDataset = "$dataset"
config.Data.inputDBS = "${DBS_URL}
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'LUMI'
config.Data.unitsPerJob = ${LUMI}

#True allows the jobs to run anywhere, regardless of where the input data is located
config.Data.ignoreLocality = True

#totalUnits only needs to be specified for GEN-SIM jobs
#config.Data.totalUnits = 200000
config.Data.outLFNDirBase = '/store/user/%s/ntuples/' % (getUsernameFromSiteDB())
config.Data.publication = True 
#config.Data.publishDataName = 'realData_FNLST_13TeV_CHNL_UNIQUE'
config.Data.outputDatasetTag =  config.General.requestName + "${productionTAG}"
config.Data.lumiMask = "$jsonFile"


#a list of the only sites at which these jobs can run
#config.Site.whitelist = ["T2_US*"]
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T2_CH_CERN'

EOF

#echo "done: $crab2File created"
crab -cfg ${crab2File} -create #|| exit 1
if [ -n "$FILELIST" ];then
    makeArguments.sh -f $FILELIST -u $UI_WORKING_DIR -n ${FILE_PER_JOB} || exit 1
fi


fi
if [ -n "${SUBMIT}" -a "`grep  DEBUG ${UI_WORKING_DIR}/log/crab.log | grep -c submit`" == "0" ];then
	echo ${UI_WORKING_DIR}
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
			mergeOutput.sh -u ${UI_WORKING_DIR} -g $file
		done
    fi
#    echo "mergeOutput.sh -u ${UI_WORKING_DIR} -n ${DATASETNAME} -r ${RUNRANGE}"
fi


done
#skim_realData_${channel[$q]}_${datasets[$q]}_${suffix[$q]}.py"



