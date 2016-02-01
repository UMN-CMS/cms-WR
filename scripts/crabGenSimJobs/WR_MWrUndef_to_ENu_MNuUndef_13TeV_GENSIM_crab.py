from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'WR_MWrUndef_to_ENu_MNuUndef_13TeV_GENSIM'
config.General.workArea = 'crab_project_WR_MWrUndef_to_ENu_MNuUndef_13TeV_GENSIM'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'WR_MWrUndef_ToENu_MNuUndef_13TeV_GEN_SIM.py'

#unitsPerJob = number of events per job
#totalUnits = total number of events to simulate
config.Data.primaryDataset = 'WR_MWrUndef_ToENu_MNuUndef_ToEEQQ_GENSIM'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 250
config.Data.totalUnits = 20000 
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.publishDataName = 'WRToENuToEEQQ_MWrUndef_MNuUndef_13TeV_GEN-SIM'

config.Site.storageSite = 'T3_US_FNALLPC'
