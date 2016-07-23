import os, sys, imp, re
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

CMSSW_VERSION=os.getenv("CMSSW_VERSION")
CMSSW_BASE=os.getenv("CMSSW_BASE")

pathPrefix=CMSSW_BASE+'/src/ExoAnalysis/cmsWR/'

config.General.requestName = 'analyze_SMPL_MMAASS_NU_MASSNU'
config.General.workArea = 'crab_project_analyze_SMPL_MMAASS_NU_MASSNU'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['output=analyzed_SMPL_MMAASS_NU_MASSNU_NUM.root']
config.JobType.psetName = pathPrefix + 'test/checkWRDecay_crabSafe_cfg.py'

#config.Data.inputDataset = 'INPTDATA'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

#use userInputFiles and primaryDataset when running over a list of input files which have not been published, but are available through xrootd
config.Data.userInputFiles = ['root://cmseos.fnal.gov//store/user/skalafut/WR/13TeV/WRSignal_March2016/WR_MWR_MMAASS_ToENu_MNu_MASSNU_GEN_13TeV_NUM.root']
config.Data.outputPrimaryDataset = 'SMPL_MMAASS_NU_MASSNU_13TeV_25ns'
#True allows the jobs to run anywhere, regardless of where the input data is located
config.Data.ignoreLocality = False

#totalUnits only needs to be specified for GEN-SIM jobs
#config.Data.totalUnits = 200000
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'analyzed_13TeV_SMPL_MMAASS_NU_MASSNU'

#a list of the only sites at which these jobs can run
#config.Site.whitelist = ["T2_US*"]
config.Site.storageSite = 'T3_US_FNALLPC'
