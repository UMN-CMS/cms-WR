import FWCore.ParameterSet.Config as cms

#import sys
#isMC=sys.modules['__main__'].isMC
#isMCsignal=sys.modules['__main__'].isMCsignal
#process = sys.modules['__main__'].process

isMC=False
isMCsignal=False
Training=False
isRun2010LoLumi=True

isData=not isMC

## Low and high lumi data selection is controlled by the JSON-derived cfi's imported
## below. For run 2010, the low lumi data is that for which the HLT_Mu9 trigger path
## was active and unprescaled, (uncertified) run range 133446 - 147116. Certification
## restricts this run range further.
##
isRun2010HiLumi=not isRun2010LoLumi

process = cms.Process("PAT");

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
# source
#process.source = cms.Source("PoolSource",
#    fileNames=cms.untracked.vstring('file:/local/cms/user/dahmes/wr2010/run2010B/oct7/mySkim/mySkim-oct7_005001.root')
#    fileNames=cms.untracked.vstring( "file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu100/HeavyNuGenHLT_WR1000_nuRmu100_1-reco-pool.root" )
#)
process.load('HeavyNu.AnalysisModules.in_cff')

if isData:
    if isRun2010LoLumi:
        print "===========> Flag is SET for LOW luminosity data <============"
        from HeavyNu.AnalysisModules.run2010loLumiRunList_cfi import lumisToProcess
    else:
        print "===========> Flag is SET for HIGH luminosity data <============"
        from HeavyNu.AnalysisModules.run2010hiLumiRunList_cfi import lumisToProcess
    
    process.source.lumisToProcess = lumisToProcess

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
if (isMC):
    print "=================> MC flag is SET <===================="
    process.GlobalTag.globaltag = cms.string('START38_V14::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_38X_V15::All')

process.load("Configuration.StandardSequences.MagneticField_cff")

################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

## pat sequences to be loaded:
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## ---
## Define the path
## ---
process.p = cms.Path(
  process.patDefaultSequence
)

########################################
# Output module - has to be defined before PAT python tools will work
########################################

if isData:
    process.out = cms.OutputModule("PoolOutputModule",
                                   fileName = cms.untracked.string('heavynu_candevents.root'),
    # save only events passing the full path
                                   SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('p')),
                                   outputCommands = cms.untracked.vstring("keep *")
                                   )
    process.outpath  = cms.EndPath(process.out)
    from PhysicsTools.PatAlgos.tools.coreTools import *
    removeMCMatching(process, ['All'])

########################################
# PAT Trigger matching
########################################
# imported directly from PhysicsTools/PatExamples/test/analyzePatTrigger_onTheFly_cfg.py
#
# PAT trigger
process.muonTriggerMatchHLTMuons = cms.EDProducer (
    "PATTriggerMatcherDRDPtLessByR",
    src            = cms.InputTag( 'selectedPatMuons' ),
    matched        = cms.InputTag( 'patTrigger' ),
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( 'TriggerMuon' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring(),
    collectionTags = cms.vstring( '*' ),
    maxDPtRel      = cms.double( 0.5 ),
    maxDeltaR      = cms.double( 0.2 ),
    maxDeltaEta    = cms.double( 0.2 ), # no effect here
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True )
)

### ============
### Python tools
### ============
### Attention: order matters!

## --
## Switch to selected PAT objects in the main work flow
## --
from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
removeCleaning( process, isData )

## --
## Switch on PAT trigger - but only for data!
## --
from PhysicsTools.PatAlgos.tools.trigTools import *
if isData:
    switchOnTriggerMatching( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ] )
    removeCleaningFromTriggerMatching( process )
    if isRun2010LoLumi: process.muonTriggerMatchHLTMuons.pathNames = cms.vstring('HLT_Mu9')
    else:               process.muonTriggerMatchHLTMuons.pathNames = cms.vstring('HLT_Mu15_v1')

##########################################
## Add analysis
##########################################

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("anal.root"),
)

process.hNu = cms.EDFilter(
    "HeavyNu",
    trigMatchPset = cms.PSet(
    trigEventTag = cms.InputTag( "" ),
    muonMatch    = cms.string( '' )
    ),
    DoLog        = cms.bool( True ),
    muonTag      = cms.InputTag( 'selectedPatMuons' ),
    jetTag       = cms.InputTag( 'selectedPatJets' ),
    electronTag  = cms.InputTag( 'selectedPatElectrons' ),
    BtagName     = cms.string('jetProbabilityBJetTags'),
    minMu1pt     = cms.double(60.),
    minMu2pt     = cms.double(15.),
    minJetPt     = cms.double(40),
    maxMu1AbsEta = cms.double(2.1),
    maxJetAbsEta = cms.double(2.5),
    minMuMuMass  = cms.double(200),
    min4objMass  = cms.double(520),
    minMuonJetdR = cms.double(0.3), # (0.5)
    muonTrackIsoLimitGeV = cms.double(10.0),

    ZmassWinMinGeV= cms.double(84.),
    ZmassWinMaxGeV= cms.double(98.),

    isSignal     = cms.bool(False),
    mNuRnormalization = cms.double(1000.0)
    )

if isData:
    # turn on trigger match requirement
    process.hNu.trigMatchPset.trigEventTag=cms.InputTag("patTriggerEvent")
    process.hNu.trigMatchPset.muonMatch=cms.string('muonTriggerMatchHLTMuons')

if Training:
    process.hNu.trainingFileName=cms.untracked.string("changeme_nntraining.txt")
    
if isMCsignal:
    process.hNu.isSignal = cms.bool(True)

process.p += process.hNu
