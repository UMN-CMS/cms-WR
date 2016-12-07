#!/bin/bash
#use this file with a command like
#bsub -R "pool>2000" -q 1nh -J runAnalysisTTMCPart_NNN < runAnalysisTTMC_NNN.sh
CMSSW_PROJECT_SRC=CMSSW_7_4_15_patch1/src/ExoAnalysis/cmsWR
UP=/afs/cern.ch/work/s/skalafut/public/WR_starting2015/wrDevelopment
export X509_USER_PROXY=/afs/cern.ch/user/s/skalafut/x509up_u38430

cd $UP/$CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
eval "cmsRun test/runAnalysis_cfg_NNN.py test=2"

