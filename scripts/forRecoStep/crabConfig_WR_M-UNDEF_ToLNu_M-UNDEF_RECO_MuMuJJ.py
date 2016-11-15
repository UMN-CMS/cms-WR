# CRAB3 config template for cmsWR

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_("General")
config.General.requestName = 'WR_M-UNDEF_ToLNu_M-UNDEF_2016-MCproductionRECO_MuMuJJ' 
config.General.workArea = 'crab_projects' 
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'WR_M-UNDEF_ToLNu_M-UNDEF_RECO_MuMuJJ.py'  

config.section_("Data")
config.Data.inputDataset = 'datasetFromDBS'
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
config.Data.outputDatasetTag = 'WR-MUNDEF_ToLNu-MUNDEF_ToMuMuJJ_RECO_13TeV-2016' 
config.Data.outLFNDirBase = '/store/user/gnegro/cmsWR/' 
config.Data.publication = True
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
config.Data.totalUnits = -1 

config.section_("Site")
config.Site.storageSite = "T2_FR_GRIF_IRFU"
config.Site.blacklist = ["T3_US_UCR", "T3_US_UMiss"] 
