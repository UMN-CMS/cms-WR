import os, sys, imp, re
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

CMSSW_VERSION=os.getenv("CMSSW_VERSION")
CMSSW_BASE=os.getenv("CMSSW_BASE")

pathPrefix=CMSSW_BASE+'/src/ExoAnalysis/cmsWR/'

config.General.requestName = 'analyze_SMPL_MMAASS_NU_NUMASS'
config.General.workArea = 'crab_project_analyze_SMPL_MMAASS_NU_NUMASS'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['output=analyzed_SMPL_MMAASS_NU_NUMASS_NUM.root']
config.JobType.psetName = pathPrefix + 'test/checkWRDecay_crabSafe_cfg.py'

#config.Data.inputDataset = 'INPTDATA'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

#use userInputFiles and primaryDataset when running over a list of input files which have not been published, but are available through xrootd
#the number of input files listed in the txt file in the line below is determined by config.Data.unitsPerJob
config.Data.userInputFiles = open(pathPrefix + "scripts/genAndFullOfflineAnalysisJobs/genWrInputFiles.txt").readlines()
config.Data.primaryDataset = 'SMPL_MMAASS_NU_NUMASS_13TeV_25ns'
#True allows the jobs to run anywhere, regardless of where the input data is located
config.Data.ignoreLocality = True

#totalUnits only needs to be specified for GEN-SIM jobs
#config.Data.totalUnits = 200000
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
#config.Data.publishDataName = 'SMPL_13TeV_analyzed_CHNL'
config.Data.outputDatasetTag = 'analyzed_13TeV_SMPL_MMAASS_NU_NUMASS'

#a list of the only sites at which these jobs can run
config.Site.whitelist = ["T2_US*"]
config.Site.storageSite = 'T3_US_FNALLPC'
