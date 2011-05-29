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
    fileNames=cms.untracked.vstring('/store/data/Run2011A/SingleMu/AOD/May10ReReco-v1/0000/0059CA13-877D-E011-BBC1-1CC1DE1CDD02.root')
)
#process.load("HeavyNu.MuFilter.in_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

## global tags: note that V16 used for Mar31, V17 used for apr18
#process.GlobalTag.globaltag = cms.string('GR_R_39X_V5::All')
#process.GlobalTag.globaltag = cms.string('GR_P_V17::All')
process.GlobalTag.globaltag = cms.string('GR_R_42_V14::All')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')


################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
 ################################################################################################

## pat sequences to be loaded:
#process.load("PhysicsTools.PFCandProducer.PF2PAT_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.patCandidateSummary.candidates = cms.VInputTag( cms.InputTag("patMuons"),cms.InputTag("patJets") )
process.patCandidates      = cms.Sequence( process.makePatMuons + process.makePatJets + process.patCandidateSummary )
process.patDefaultSequence = cms.Sequence( process.patCandidates )

process.muFilter = cms.EDFilter("MuFilter",
    minPt    = cms.double( 15.0 ),
    keepJets = cms.bool( False ),
    #--- Included, but we ignore this parameter ---#
    HLTpaths = cms.vstring( 'HLT_Mu9','HLT_Mu11' )
)

process.mySkimPath = cms.Path(
    process.patDefaultSequence *
    # process.hNu *
    process.muFilter 
)

process.out = cms.OutputModule( "PoolOutputModule",
   fileName = cms.untracked.string("poolout.root"),
   maxSize = cms.untracked.int32(3000000),
   SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('mySkimPath'),
   ),
   outputCommands = cms.untracked.vstring("keep *")
)

process.outpath = cms.EndPath(process.out)

from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData(process, ['All'])

#
# kluge for 42X, because L2L3Residual corrections not yet available.
#
if 'L2L3Residual' in process.patJetCorrFactors.levels:
   process.patJetCorrFactors.levels.remove('L2L3Residual')
