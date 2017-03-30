#!/bin/bash
#use this file with a command like
#bsub -R "pool>2000" -q 1nh -J jobName_NNN < jobRunScript_NNN.sh
CMSSW_PROJECT_SRC=CMSSW_8_0_8_patch1/src/ExoAnalysis/cmsWR
UP=/afs/cern.ch/work/s/skalafut/public/WR_starting2015/wrDevelopment/effWithOffDiagonalSamples
export X509_USER_PROXY=/afs/cern.ch/user/s/skalafut/x509up_u38430

cd $UP/$CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
eval "cmsRun test/runAnalysis_cfg_TAGNAME_NNN.py test=2"

