#!/bin/bash
#use this file with a command like
#bsub -R "rusage[mem=8000]" -q 1nd -J runAnalysisTTMC < runAnalysisTTMC.sh
CMSSW_PROJECT_SRC=CMSSW_7_4_15_patch1/src/ExoAnalysis/cmsWR
UP=/afs/cern.ch/work/s/skalafut/public/WR_starting2015/wrDevelopment
export X509_USER_PROXY=/afs/cern.ch/user/s/skalafut/x509up_u38430

cd $UP/$CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
eval "cmsRun test/runAnalysis_cfg.py test=2 output=reprocessTTJetsV2MCMinitrees.root >& TTMCOutput.txt &"

