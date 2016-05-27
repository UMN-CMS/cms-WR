#!/bin/bash
#use this file with a command like
#bsub -R "rusage[mem=8000]" -q 8nh -J wrEE < runWrSignalEEJob.sh
#to request 8 GB ram
CMSSW_PROJECT_SRC=CMSSW_7_4_15_patch1/src/ExoAnalysis/cmsWR
UP=/afs/cern.ch/user/s/skalafut/WR
finalDir=/afs/cern.ch/work/s/skalafut/public/WR_starting2015/forPeterRootFiles
export X509_USER_PROXY=/afs/cern.ch/user/s/skalafut/x509up_u38430

cd $UP/$CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
eval 'make'
eval "./bin/analysis -m signal -c EE --toys 200 --nStatToys 200"
eval "mv ./selected_tree_WRtoEEJJ_*signal_eeEE.root $finalDir/."

