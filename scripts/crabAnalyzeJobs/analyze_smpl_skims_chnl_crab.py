from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'analyze_SMPL_skim_low_mass_region_CHNL'
config.General.workArea = 'crab_project_analyze_SMPL_skim_low_mass_region_CHNL'
config.General.transferOutputs = True
config.General.transferLogs = True


config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/uscms/home/skalafut/nobackup/WR_starting2015/mostUpToDateCode/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/test/quickRecoKinematicsAnalyzer_CHNL_cfg.py'
#config.JobType.maxMemoryMB = 2500 should not need this option for skims

config.Data.inputDataset = 'INPTDATA'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#True allows the jobs to run anywhere, regardless of where the input data is located
#config.Data.ignoreLocality = True

#totalUnits only needs to be specified for GEN-SIM jobs
#config.Data.totalUnits = 200000
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True 
config.Data.publishDataName = 'SMPL_13TeV_25ns_analyzed_CHNL_skims_low_mass_region'

#config.Site.whitelist = ["T2_US*"] a list of the only sites at which these jobs can run
config.Site.storageSite = 'T3_US_FNALLPC'
