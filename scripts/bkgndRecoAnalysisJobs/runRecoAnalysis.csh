#!/bin/csh

#Setup environment, move and untar needed files and directories, compile everything
setenv HOME /eos/uscms/store/user/skalafut/WR/13TeV/bkgnds/BBKKGGNNDD
setenv PATH /bin:/usr/bin:/usr/local/bin:/usr/krb5/bin:/usr/afsws/bin:/usr/krb5/bin/aklog

#source /uscmst1/prod/sw/cms/cshrc prod
source /cvmfs/cms.cern.ch/cmsset_default.csh prod
setenv SCRAM_ARCH slc6_amd64_gcc481

cmsrel CMSSW_7_4_0_pre9

mv forBkgndAnalysis.tar bkgndRecoAnalysis_BBKKGGNNDD_NUM.py CMSSW_7_4_0_pre9/src/.

cd CMSSW_7_4_0_pre9/src/
tar -xvf forBkgndAnalysis.tar
rm forBkgndAnalysis.tar

cmsenv
eval `scramv1 runtime -csh`
eval 'scram b -j 8'

cmsRun bkgndRecoAnalysis_BBKKGGNNDD_NUM.py

mv analyzed_BBKKGGNNDD_bkgndElectronChannel_NUM.root $HOME/.
 
