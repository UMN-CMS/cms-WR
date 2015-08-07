from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'FNLST_skim_low_mass_and_signal_regions_MINIAOD_atFNALLPC'
config.General.workArea = 'crab_project_FNLST_skim_low_mass_and_signal_regions_MINIAOD_atFNALLPC'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/uscms/home/skalafut/nobackup/WR_starting2015/mostUpToDateCode/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/test/skim_bkgndMC_CHNL_signal_and_low_mass_regions_cfg.py'
#config.JobType.maxMemoryMB = 2500 should not need this option for skims

config.Data.inputDataset = 'INPTDT'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
#config.Data.ignoreLocality = True allows the jobs to run anywhere, regardless of where the input data is located

#totalUnits only needs to be specified for GEN-SIM jobs
#config.Data.totalUnits = 200000
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.publishDataName = 'FNLST_13TeV_skim_low_mass_and_signal_regions_atFNALLPC'

#a list of the only sites at which these jobs can run
#config.Site.whitelist = ["T2_US*"]
#config.Site.storageSite = 'T3_US_Cornell'
config.Site.storageSite = 'T3_US_FNALLPC'

