from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'FNLST_July2015_CHNL_signalAndLowMassRegionSkim'
config.General.workArea = 'crab_project_FNLST_July2015_CHNL_signalAndLowMassRegionSkim'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/uscms/home/skalafut/nobackup/WR_starting2015/mostUpToDateCode/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/test/skim_realData_CHNL_signal_and_low_mass_regions_cfg.py'
#config.JobType.maxMemoryMB = 2500 should not need this option for skims

config.Data.inputDataset = '/FNLST/Run2015B-PromptReco-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/uscms/home/skalafut/nobackup/WR_starting2015/crabDir/realData/Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt'
config.Data.unitsPerJob = 30000
#True allows the jobs to run anywhere, regardless of where the input data is located
#config.Data.ignoreLocality = True

#totalUnits only needs to be specified for GEN-SIM jobs
#config.Data.totalUnits = 200000
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True 
config.Data.publishDataName = 'realData_FNLST_13TeV_50ns_CHNL_signalAndLowMassRegionSkim'

#config.Site.whitelist = ["T2_US*"] a list of the only sites at which these jobs can run
config.Site.storageSite = 'T3_US_Cornell'
