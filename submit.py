import os

string = '''
from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName   = 'runAnalysis_80X_%(name)s'
config.General.transferLogs = True
config.General.workArea = 'crab'
config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'test/runAnalysis_cfg.py'
config.JobType.pyCfgParams = ['isMC=%(MC)s','datasetTag=%(dname)s']
config.JobType.inputFiles = ['MCPileup.root','DataPileup.root','python/MuonHighPt_Z_RunCD_Reco74X_Dec17.pkl']

config.section_("Data")
config.Data.inputDataset = '%(dataset)s'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10

config.Data.outLFNDirBase = '/store/user/%(user)s/runAnalysis_80X_%(name)s'


config.Data.lumiMask = '%(json)s'

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'

'''

user = 'jchavesb'

datasets = []
datasetTags = []

f = open('configs/datasets_80X.dat')
for line in f:
    if '#' not in line:
        datasetTags.append(line.split('\t')[0])
        datasets.append(line.split('\t')[1])

json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON.txt'

for d,dt in zip(datasets,datasetTags):    
    vd = locals()
    vd['user'] = user
    vd['dataset'] = d
    vd['name'] = dt
    vd['MC'] = '1'
    vd['json'] = ''
    if 'Run2016' in d:
        vd['MC'] = '0'
        vd['json'] = json
    vd['dname'] = dt
    open('crabConfig_tmp.py','wt').write(string % vd)
    os.system('crab submit -c crabConfig_tmp.py')
