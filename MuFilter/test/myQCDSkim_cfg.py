import FWCore.ParameterSet.Config as cms

process = cms.Process("QCDskim")

#--- Filter will collect electron QCD or muon QCD, not both ---#
doMuons     = True
doElectrons = not doMuons

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
# source
process.source = cms.Source("PoolSource",      
    fileNames=cms.untracked.vstring('file:input.root')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# process.GlobalTag.globaltag = cms.string('GR_R_50_V13::All')
process.GlobalTag.globaltag = cms.string('FT_53_V6_AN1::All')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')

### Output needs to be created before working with PAT objects ###
process.out = cms.OutputModule( "PoolOutputModule",
   fileName = cms.untracked.string("poolout.root"),
   maxSize = cms.untracked.int32(3000000),
   SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('p'),
   ),
   outputCommands = cms.untracked.vstring("keep *","drop *_*_*_QCDskim","keep edmTriggerResults_TriggerResults__QCDskim")
)

###########################################################################
###  P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s  ###
###########################################################################

# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")
# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences. 
# from PhysicsTools.PatAlgos.tools.pfTools import *
# postfix = "PFlow"
# jetAlgo = "AK5"
# usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=False, postfix=postfix) 

process.patCandidateSummary.candidates = cms.VInputTag( cms.InputTag("patMuons"), cms.InputTag("patElectrons") )
process.patCandidates      = cms.Sequence( process.makePatMuons + process.makePatElectrons + process.patCandidateSummary )
process.patDefaultSequence = cms.Sequence( process.patCandidates )

process.qcdFilter = cms.EDFilter("QCDFilter",
    muonTag     = cms.InputTag("patMuons"),
    electronTag = cms.InputTag("patElectrons"),

    collectMuons     = cms.bool(doMuons),
    collectElectrons = cms.bool(doElectrons),

    minMuonPt     = cms.double( 35.0 ),
    minElectronPt = cms.double( 35.0 ),

    lowPtPrescale = cms.int32( 5 ),
    lowPtCutoff   = cms.double( 60.0 ),
)

process.hltMuonFilter = cms.EDFilter("TriggerResultsFilter",
     hltResults   = cms.InputTag('TriggerResults','','HLT'), # HLT results   - set to empty to ignore HLT
     l1tResults   = cms.InputTag(''),                        # L1 GT results - set to empty to ignore L1
     l1tIgnoreMask         = cms.bool(False),   # use L1 mask
     l1techIgnorePrescales = cms.bool(False),   # read L1 technical bits from PSB#9, bypassing the prescales
     daqPartitions         = cms.uint32(0x01),  # used by the definition of the L1 mask
     throw                 = cms.bool(False),   # if HLT path not in the table, crash/ignore according to true/false
     triggerConditions     = cms.vstring(
           'HLT_Mu12_v16', ## Needed only up to run 193834, when HLT_Mu40 was activated and unprescaled 
           # 'HLT_Mu12_v17', ## Not needed as HLT_Mu40 activated and unprescaled by this time
           # 'HLT_Mu12_v18', ## Not needed as HLT_Mu40 activated and unprescaled by this time
           'HLT_Mu40_v12',
           'HLT_Mu40_v13',
           'HLT_Mu40_v14',
           'HLT_Mu40_eta2p1_v9',
           'HLT_Mu40_eta2p1_v10',
           'HLT_Mu40_eta2p1_v11'
     )
)

process.hltMuon12Filter = process.hltMuonFilter.clone()
process.hltMuon12Filter.triggerConditions = cms.vstring( 
    'HLT_Mu12_v16', ## Needed only up to run 193834, when HLT_Mu40 was activated and unprescaled 
    'HLT_Mu12_v17', ## Not needed as HLT_Mu40 activated and unprescaled by this time
    'HLT_Mu12_v18' ## Not needed as HLT_Mu40 activated and unprescaled by this time
)
process.hltMuon40Filter = process.hltMuonFilter.clone()
process.hltMuon40Filter.triggerConditions = cms.vstring( 
    'HLT_Mu40_v12',
    'HLT_Mu40_v13',
    'HLT_Mu40_v14'
)
process.hltMuon40eta2p1Filter = process.hltMuonFilter.clone()
process.hltMuon40eta2p1Filter.triggerConditions = cms.vstring( 
    'HLT_Mu40_eta2p1_v9',
    'HLT_Mu40_eta2p1_v10',
    'HLT_Mu40_eta2p1_v11'
)

process.hltElectronFilter = process.hltMuonFilter.clone()
process.hltElectronFilter.triggerConditions = cms.vstring( 
    'HLT_Photon30_CaloIdVL_v11',
    'HLT_Photon30_CaloIdVL_v12',
    'HLT_Photon30_CaloIdVL_v13',
    'HLT_Photon30_CaloIdVL_v14',
    'HLT_Photon50_CaloIdVL_v7',
    'HLT_Photon50_CaloIdVL_v8',
    'HLT_Photon50_CaloIdVL_v9',
    'HLT_Photon50_CaloIdVL_v10',
    'HLT_Photon75_CaloIdVL_v10', 
    'HLT_Photon75_CaloIdVL_v11', 
    'HLT_Photon75_CaloIdVL_v12', 
    'HLT_Photon75_CaloIdVL_v13', 
    'HLT_Photon90_CaloIdVL_v7',
    'HLT_Photon90_CaloIdVL_v8',
    'HLT_Photon90_CaloIdVL_v9',
    'HLT_Photon90_CaloIdVL_v10',
    'HLT_Photon135_CaloIdVL_v4',
    'HLT_Photon135_CaloIdVL_v5',
    'HLT_Photon135_CaloIdVL_v6',
    'HLT_Photon135_CaloIdVL_v7',
    'HLT_Photon150_CaloIdVL_v1',
    'HLT_Photon150_CaloIdVL_v2',
    'HLT_Photon150_CaloIdVL_v3',
    'HLT_Photon150_CaloIdVL_v4', 
    'HLT_Photon160_CaloIdVL_v1', ## Unnecessary as HLT_Photon150 unprescaled everywhere
    'HLT_Photon160_CaloIdVL_v2', ## Unnecessary as HLT_Photon150 unprescaled everywhere
    'HLT_Photon160_CaloIdVL_v3', ## Unnecessary as HLT_Photon150 unprescaled everywhere
    'HLT_Photon160_CaloIdVL_v4'  ## Unnecessary as HLT_Photon150 unprescaled everywhere
)

# The following are special filters to verify HLT paths (i.e. protect against typos)
process.hltPhoton30Filter = process.hltElectronFilter.clone()
process.hltPhoton30Filter.triggerConditions = cms.vstring( 
    'HLT_Photon30_CaloIdVL_v11',
    'HLT_Photon30_CaloIdVL_v12',
    'HLT_Photon30_CaloIdVL_v13',
    'HLT_Photon30_CaloIdVL_v14'
)
process.hltPhoton50Filter = process.hltElectronFilter.clone()
process.hltPhoton50Filter.triggerConditions = cms.vstring( 
    'HLT_Photon50_CaloIdVL_v7',
    'HLT_Photon50_CaloIdVL_v8',
    'HLT_Photon50_CaloIdVL_v9',
    'HLT_Photon50_CaloIdVL_v10'
)
process.hltPhoton75Filter = process.hltElectronFilter.clone()
process.hltPhoton75Filter.triggerConditions = cms.vstring( 
    'HLT_Photon75_CaloIdVL_v10', 
    'HLT_Photon75_CaloIdVL_v11', 
    'HLT_Photon75_CaloIdVL_v12', 
    'HLT_Photon75_CaloIdVL_v13'
)
process.hltPhoton90Filter = process.hltElectronFilter.clone()
process.hltPhoton90Filter.triggerConditions = cms.vstring( 
    'HLT_Photon90_CaloIdVL_v7',
    'HLT_Photon90_CaloIdVL_v8',
    'HLT_Photon90_CaloIdVL_v9',
    'HLT_Photon90_CaloIdVL_v10'
)
process.hltPhoton135Filter = process.hltElectronFilter.clone()
process.hltPhoton135Filter.triggerConditions = cms.vstring( 
    'HLT_Photon135_CaloIdVL_v4',
    'HLT_Photon135_CaloIdVL_v5',
    'HLT_Photon135_CaloIdVL_v6',
    'HLT_Photon135_CaloIdVL_v7'
)
process.hltPhoton150Filter = process.hltElectronFilter.clone()
process.hltPhoton150Filter.triggerConditions = cms.vstring( 
    'HLT_Photon150_CaloIdVL_v1',
    'HLT_Photon150_CaloIdVL_v2',
    'HLT_Photon150_CaloIdVL_v3',
    'HLT_Photon150_CaloIdVL_v4'
)

if doMuons:
    process.p = cms.Path( process.hltMuonFilter * process.patDefaultSequence * process.qcdFilter )
    process.pMuon12       = cms.Path( process.hltMuon12Filter ) 
    process.pMuon40       = cms.Path( process.hltMuon40Filter ) 
    process.pMuon40eta2p1 = cms.Path( process.hltMuon40eta2p1Filter ) 
else:
    process.p = cms.Path( process.hltElectronFilter * process.patDefaultSequence * process.qcdFilter )
    process.pPhoton30  = cms.Path( process.hltPhoton30Filter ) 
    process.pPhoton50  = cms.Path( process.hltPhoton50Filter ) 
    process.pPhoton75  = cms.Path( process.hltPhoton75Filter ) 
    process.pPhoton90  = cms.Path( process.hltPhoton90Filter ) 
    process.pPhoton135 = cms.Path( process.hltPhoton135Filter ) 
    process.pPhoton150 = cms.Path( process.hltPhoton150Filter ) 

process.outpath = cms.EndPath(process.out)

from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData(process, ['All'])

#
# kluge for 42X, because L2L3Residual corrections not yet available.
#
# if 'L2L3Residual' in process.patJetCorrFactors.levels:
#    process.patJetCorrFactors.levels.remove('L2L3Residual')
