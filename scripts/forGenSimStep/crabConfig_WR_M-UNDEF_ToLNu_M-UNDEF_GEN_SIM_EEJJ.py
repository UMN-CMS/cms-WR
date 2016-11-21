# CRAB3 config template for cmsWR

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_("General")
config.General.requestName = 'WR_M-UNDEF_ToLNu_M-UNDEF_2016-MCproduction_EEJJ' 
config.General.workArea = 'crab_projects'  
config.General.transferOutputs = True
config.General.transferLogs = True  #False

config.section_("JobType")
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'WR_M-UNDEF_ToLNu_M-UNDEF_GEN_SIM_EEJJ.py'  

config.section_("Data")
config.Data.outputPrimaryDataset = 'WR-ToLNu-ToEEJJ_GEN_SIM_13TeV-2016' 
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 250 
NJOBS = 200
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/gnegro/cmsWR/' 
config.Data.publication = True
config.Data.outputDatasetTag = 'WR-MUNDEF_ToLNu-MUNDEF_ToEEJJ_GEN_SIM_13TeV-2016' 
# config.Data.publishDBS = 'phys03'  

config.section_("Site")
config.Site.storageSite = "T2_FR_GRIF_IRFU"
config.Site.blacklist = ["T3_US_UCR", "T2_US_Vanderbilt", "T2_EE_Estonia"] 
