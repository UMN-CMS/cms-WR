from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'FNLST_CHNL_UNIQUE'
config.General.workArea = 'crab_project_FNLST_CHNL_UNIQUE'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/uscms/home/skalafut/nobackup/WR_starting2015/mostUpToDateCode/Run2015D/CMSSW_7_4_12_patch2/src/ExoAnalysis/cmsWR/test/skim_realData_CHNL_signal_and_low_mass_regions_cfg.py'
#config.JobType.psetName = '/uscms/home/skalafut/nobackup/WR_starting2015/mostUpToDateCode/Run2015C/CMSSW_7_4_10_patch1/src/ExoAnalysis/cmsWR/test/skim_realData_CHNL_signal_and_low_mass_regions_cfg.py'
#config.JobType.maxMemoryMB = 2500 should not need this option for skims

config.Data.inputDataset = '/FNLST/TAG/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'LUMI'
config.Data.unitsPerJob = 30 

#True allows the jobs to run anywhere, regardless of where the input data is located
config.Data.ignoreLocality = True

#totalUnits only needs to be specified for GEN-SIM jobs
#config.Data.totalUnits = 200000
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True 
#config.Data.publishDataName = 'realData_FNLST_13TeV_CHNL_UNIQUE'
config.Data.outputDatasetTag = 'realData_FNLST_13TeV_CHNL_UNIQUE'

#a list of the only sites at which these jobs can run
config.Site.whitelist = ["T2_US*"]
config.Site.storageSite = 'T3_US_FNALLPC'
