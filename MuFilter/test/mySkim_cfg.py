import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
# source
process.source = cms.Source("PoolSource",      
    fileNames=cms.untracked.vstring(  )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
process.GlobalTag.globaltag = cms.string('GR_R_36X_V12::All')
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
# process.patElectrons.addGenMatch  = cms.bool(False)
process.patMuons.addGenMatch      = cms.bool(False)
# process.patTaus.addGenMatch       = cms.bool(False)
# process.patTaus.addGenJetMatch    = cms.bool(False)
# process.patPhotons.addGenMatch    = cms.bool(False)
# process.patJets.addGenPartonMatch = cms.bool(False)
# process.patJets.addGenJetMatch    = cms.bool(False)
# process.patMETs.addGenMET         = cms.bool(False)
# process.patElectrons.embedGenMatch  = cms.bool(False)
process.patMuons.embedGenMatch      = cms.bool(False)
# process.patTaus.embedGenMatch       = cms.bool(False)
# process.patTaus.embedGenJetMatch    = cms.bool(False)
# process.patPhotons.embedGenMatch    = cms.bool(False)
# process.patJets.embedGenPartonMatch = cms.bool(False)
# process.patJets.embedGenJetMatch    = cms.bool(False)
process.makePatMuons  = cms.Sequence( process.patMuons )
process.patCandidateSummary.candidates = cms.VInputTag( cms.InputTag("patMuons") )

process.patCandidates      = cms.Sequence( process.makePatMuons + process.patCandidateSummary )
process.patDefaultSequence = cms.Sequence( process.patCandidates )

# process.makePatElectrons = cms.Sequence( process.patElectrons )
# process.makePatTaus      = cms.Sequence( process.patTaus )
# process.makePatPhotons   = cms.Sequence( process.patPhotons )
# process.makePatJets      = cms.Sequence( process.patJets )

# THIS IS PROD
process.TFileService = cms.Service("TFileService",
       fileName = cms.string(''),
)

process.muFilter = cms.EDFilter("MuFilter",
    minPt = cms.double( 15.0 ),
    HLTpaths = cms.vstring( 'HLT_Mu9','HLT_Mu11' )
)

process.hNu = cms.EDAnalyzer("HeavyNu",
       DoLog = cms.bool( True )
)

process.mySkimPath = cms.Path(
    process.patDefaultSequence *
    # process.hNu *
    process.muFilter 
)

process.wrOutputModule = cms.OutputModule( "PoolOutputModule",
   fileName = cms.untracked.string(""),
   maxSize = cms.untracked.int32(300000),
   SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('mySkimPath'),
   )
)

process.outpath = cms.EndPath(process.wrOutputModule)
