from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'WR_MWrUndef_to_LNu_MNuUndef_13TeV_NNUUMMPU_DIGIRAWRECO'
config.General.workArea = 'crab_project_WR_MWrUndef_to_LNu_MNuUndef_13TeV_NNUUMMPU_DIGIRAWRECO'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'WR_MWrUndef_ToLNu_MNuUndef_13TeV_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_NNUUMMPU.py'
config.JobType.maxMemoryMB = 3500


config.Data.inputDataset = 'datasetFromDBS'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = nInputFiles

#this line is needed for running over files published to FNALLPC T3 site
config.Data.ignoreLocality = True

config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.publishDataName = 'WRToLNu_ToLLJJ_MWrUndef_MNuUndef_13TeV_NNUUMMPU_RAW_RECO'

config.Site.whitelist = ["T2_US*"]
config.Site.storageSite = 'publishSite'
