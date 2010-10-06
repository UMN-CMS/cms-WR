import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
# source
process.source = cms.Source("PoolSource",      
    fileNames=cms.untracked.vstring( "file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu300/HeavyNuGenHLT_WR1000_nuRmu300_1-reco-pool.root" )                 
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
process.GlobalTag.globaltag = cms.string('START36_V10::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')


################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

## pat sequences to be loaded:
#process.load("PhysicsTools.PFCandProducer.PF2PAT_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
#process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
##

# THIS IS PROD
process.TFileService = cms.Service("TFileService",
       fileName = cms.string("foo.root"),
)

process.hNu = cms.EDAnalyzer("HeavyNu",
       DoLog = cms.bool( True )
)

process.myAnalysisPath = cms.Path(
    process.patDefaultSequence *
    process.hNu 
)
