from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'analyze_SMPL'
config.General.workArea = 'crab_project_analyze_SMPL'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/uscms/home/skalafut/nobackup/WR_starting2015/mostUpToDateCode/CMSSW_7_4_12_patch4/src/ExoAnalysis/cmsWR/test/unmatched_recoElectronChannel_noCmdLineInputs_cfg.py'
#config.JobType.psetName = '/uscms/home/skalafut/nobackup/WR_starting2015/mostUpToDateCode/CMSSW_7_4_12_patch4/src/ExoAnalysis/cmsWR/test/checkWRDecay_crabSafe_cfg.py'
#config.JobType.maxMemoryMB = 2500 should not need this option for skims

config.Data.inputDataset = 'INPTDATA'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 4500

#use userInputFiles and primaryDataset when running over a list of input files which have not been published, but are available through xrootd
#config.Data.userInputFiles = open('').readlines()
#config.Data.primaryDataset = 'SMPL_13TeV_25ns_CHNL_skim_low_dilepton_mass_region'
#True allows the jobs to run anywhere, regardless of where the input data is located
#config.Data.ignoreLocality = True

#totalUnits only needs to be specified for GEN-SIM jobs
#config.Data.totalUnits = 200000
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
#config.Data.publishDataName = 'SMPL_13TeV_analyzed_CHNL'
config.Data.outputDatasetTag = 'SMPL_13TeV_analyzed_CHNL'

#a list of the only sites at which these jobs can run
#config.Site.whitelist = ["T2_US*"]
config.Site.storageSite = 'T3_US_FNALLPC'
