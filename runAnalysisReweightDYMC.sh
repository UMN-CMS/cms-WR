#!/bin/bash
#use this file with a command like
#bsub -R "pool>2000" -q 8nh -J uniqueJobName_WGTNUM < runAnalysisReweightDYMC_MODE_CHNL_WGTNUM.sh
CMSSW_PROJECT_SRC=CMSSW_7_4_15_patch1/src/ExoAnalysis/cmsWR
UP=/afs/cern.ch/work/s/skalafut/public/WR_starting2015/wrDevelopment
#export X509_USER_PROXY=/afs/cern.ch/user/s/skalafut/x509up_u38430

cd $UP/$CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
eval "make"
#output file size is too large with tagandprobe eval "./bin/analysis -m MODE -c CHNL --isTagAndProbe true --pdfWeight WGTNUM"
#eval "./bin/analysis -m MODE -c CHNL --ignoreDyScaleFactors false --isLowDiLepton true --pdfWeight WGTNUM"
eval "./bin/analysis -m MODE -c CHNL --ignoreDyScaleFactors false --pdfWeight WGTNUM"

