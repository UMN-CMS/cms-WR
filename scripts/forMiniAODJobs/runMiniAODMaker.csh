#!/bin/csh

#Setup environment, move and untar needed files and directories, compile everything
setenv HOME /eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_ABCDE
setenv PATH /bin:/usr/bin:/usr/local/bin:/usr/krb5/bin:/usr/afsws/bin:/usr/krb5/bin/aklog

#source /uscmst1/prod/sw/cms/cshrc prod
source /cvmfs/cms.cern.ch/cmsset_default.csh prod
setenv SCRAM_ARCH slc6_amd64_gcc481

cmsrel CMSSW_7_4_0_pre9

mv PU.tar genSimFile_NUM.root reco_and_miniAOD_NUM.py CMSSW_7_4_0_pre9/src/.

cd CMSSW_7_4_0_pre9/src/
tar -xvf PU.tar
rm PU.tar

cmsenv
eval `scramv1 runtime -csh`
eval 'scram b -j 8'

cmsRun reco_and_miniAOD_NUM.py

mv ABCDE_miniAODFile_NUM.root $HOME/.
 
