import os,sys
import FWCore.ParameterSet.Config as cms

process = cms.Process('ANA')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('ExoAnalysis.cmsWR.treeMaker_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('file:/afs/cern.ch/user/j/jchavesb/work/WR_skims/DYJets/skim_signal_1.root'),
                            )

os.system('das_client.py --query="file dataset=/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/jchavesb-crab_WR_DY_Jets_Skim_50ns_Sideband-8a626735db4a89bf68c59f2da6d5e98a/USER instance=prod/phys03" --limit=0 > pappy.txt')

fn = 'pappy.txt'

sec_file = open(fn, 'r')

mysecfilelist = []
for line in sec_file:
    # add as many files as you wish this way
    if 'skim' in line:
        #mysecfilelist.append('file:'+line.strip())
        mysecfilelist.append(line.strip())
#process.source.fileNames = mysecfilelist
outfn = 'skim_ttree.root'
process.TFileService = cms.Service('TFileService', fileName = cms.string(outfn))

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)


from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.triggerFilter = hltHighLevel.clone()
process.triggerFilter.HLTPaths = ['HLT_Mu45_eta2p1_v*','HLT_Mu50_v*','HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*']
process.triggerFilter.andOr = True # = OR
#for name, path in process.paths.items():
 #   if not name.startswith('eventCleaning'):
  #      path.insert(0, process.triggerFilter)
process.ptrig = cms.Path(process.triggerFilter)

process.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('ptrig'))

process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Summer15_50nsV2_MC_AK4PFchs'),
            # tag    = cms.string('JetCorrectorParametersCollection_Summer12_V3_MC_AK5PF'),
            label  = cms.untracked.string('AK4PFchs')
            ),
      ## here you add as many jet types as you need
      ## note that the tag name is specific for the particular sqlite file 
      ), 
      connect = cms.string('sqlite:Summer15_50nsV2_MC.db')
)
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')


from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
  src = cms.InputTag("slimmedJets"),
  levels = cms.vstring('L1FastJet', 
        'L2Relative', 
        'L3Absolute'),
  payload = cms.string('AK4PFchs') ) # Make sure to choose the appropriate levels and payload here!

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = patJetsUpdated.clone(
  jetSource = cms.InputTag("slimmedJets"),
  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
  )

#process.p += cms.Sequence( process.patJetCorrFactorsReapplyJEC + process. patJetsReapplyJEC )


#process.p = cms.Path(process.patJetCorrFactorsReapplyJEC + process. patJetsReapplyJEC + process.wRsubleadingElectron + process.wRlooseJet + process.triggerFilter  + process.MakeTTree_Electrons)
process.p = cms.Path(process.patJetCorrFactorsReapplyJEC + process. patJetsReapplyJEC + process.wRtunePMuons + process.wRsubleadingMuon + process.wRlooseJet + process.triggerFilter  + process.MakeTTree_Muons)

