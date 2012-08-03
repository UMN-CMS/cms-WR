import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

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

process.GlobalTag.globaltag = cms.string('GR_R_50_V13::All')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')

### Output needs to be created before working with PAT objects ###
process.out = cms.OutputModule( "PoolOutputModule",
   fileName = cms.untracked.string("poolout.root"),
   maxSize = cms.untracked.int32(3000000),
   SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('p'),
   ),
   outputCommands = cms.untracked.vstring("keep *","drop *_*_*_SKIM","keep edmTriggerResults_TriggerResults__SKIM")
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

process.patCandidateSummary.candidates = cms.VInputTag( cms.InputTag("patMuons") )
process.patCandidates      = cms.Sequence( process.makePatMuons + process.patCandidateSummary )
process.patDefaultSequence = cms.Sequence( process.patCandidates )

#--- Generic PAT tracks modules stolen from ElectroWeakAnalysis/Skimming/python ---#
process.patAODTrackCandsUnfiltered = cms.EDProducer("ConcreteChargedCandidateProducer",
    src          = cms.InputTag("generalTracks"),
    particleType = cms.string('mu+')   # to fix mass hypothesis
)
process.patAODTrackCands = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("patAODTrackCandsUnfiltered"),
    cut = cms.string('pt > 20')
)
from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import patGenericParticles
process.allPatTracks = patGenericParticles.clone(
    src = cms.InputTag("patAODTrackCands")
)
from PhysicsTools.PatAlgos.selectionLayer1.trackSelector_cfi import *
process.patTracksPt20 = selectedPatTracks.clone(
    cut = 'pt > 30.'
)

process.muFilter = cms.EDFilter("MuFilter",
    muonTag  = cms.InputTag("patMuons"),
    trackTag = cms.InputTag("patTracksPt20"),
    ebTag    = cms.InputTag("correctedHybridSuperClusters"),
    eeTag    = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),

    minMuonPt = cms.double( 32.0 ),
    minSCEt   = cms.double( 32.0 ), 
    overlap   = cms.double( 0.05 ),
    trackPrescale = cms.int32( -1 ),
    trackOnly = cms.bool(False),
    trackMassLow = cms.double( 40.0),
    trackMassHigh = cms.double( 150.0)
)

process.patTrackSequence = cms.Sequence( 
        process.patAODTrackCandsUnfiltered *
        process.patAODTrackCands *
        process.allPatTracks *
        process.patTracksPt20
)

process.p = cms.Path( 
   # getattr(process,"patPF2PATSequence"+postfix) * 
   process.patDefaultSequence * 
   process.patTrackSequence * 
   process.muFilter
)

process.outpath = cms.EndPath(process.out)

from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData(process, ['All'])

#
# kluge for 42X, because L2L3Residual corrections not yet available.
#
# if 'L2L3Residual' in process.patJetCorrFactors.levels:
#    process.patJetCorrFactors.levels.remove('L2L3Residual')
