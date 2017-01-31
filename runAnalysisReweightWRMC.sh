#!/bin/bash
#use this file with a command like
#bsub -R "pool>1500" -q 1nh -J uniqueJobName_WGTNUM < runAnalysisReweightWRMC_MODE_CHNL_WGTNUM.sh
export X509_USER_PROXY=/afs/cern.ch/user/s/skalafut/x509up_u38430

cd WORKDIR
eval `scramv1 runtime -sh`
eval "make"
eval "./bin/analysis -m MODE -c CHNL --pdfWeight WGTNUM --signalN MMAASS"

