#!/bin/bash
#use this file with a command like
#bsub -R "pool>2000" -q 1nh -J runAnalysisDYMCPart_NNN < runAnalysisDYMC_NNN.sh
CMSSW_PROJECT_SRC=CMSSW_7_4_15_patch1/src/ExoAnalysis/cmsWR
UP=/afs/cern.ch/work/s/skalafut/public/WR_starting2015/wrDevelopment/forEleSystematic
export X509_USER_PROXY=/afs/cern.ch/user/s/skalafut/x509up_u38430
WRITEDIR='/eos/cms/store/group/phys_exotica/leptonsPlusJets/WR/skims/TAGNAME_SHv12'
eosReadingTag='root://eoscms.cern.ch/'

cd $UP/$CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
#dont save TnP branch, it consumes too much disk space
eval "cmsRun test/skims_cfg_TAGNAME_NNN.py"

#after the skim finishes, move the skim file to a pre-existing directory in the WR group space on eos
eval "xrdcp TAGNAME_skimPartNNN.root $eosReadingTag$WRITEDIR/."
rm TAGNAME_skimPartNNN.root

