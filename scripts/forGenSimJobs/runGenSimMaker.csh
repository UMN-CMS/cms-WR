#!/bin/csh

#Setup environment, move and untar needed files and directories, compile everything
setenv HOME /eos/uscms/store/user/skalafut/WR/13TeV/Signal_GEN-SIM
setenv PATH /bin:/usr/bin:/usr/local/bin:/usr/krb5/bin:/usr/afsws/bin:/usr/krb5/bin/aklog

#source /uscmst1/prod/sw/cms/cshrc prod
source /cvmfs/cms.cern.ch/cmsset_default.csh prod
setenv SCRAM_ARCH slc6_amd64_gcc481

cmsrel CMSSW_7_4_0_pre9

mv WR_M-MMAASS_ToLNu_M-MASSNU_GEN_SIM_NUM.py CMSSW_7_4_0_pre9/src/.

cd CMSSW_7_4_0_pre9/src/

cmsenv
eval `scramv1 runtime -csh`
eval 'scram b -j 8'

cmsRun WR_M-MMAASS_ToLNu_M-MASSNU_GEN_SIM_NUM.py

mv WR_M-MMAASS_ToLNu_M-MASSNU_GEN_SIM_13TeV_NUM.root $HOME/.
 
