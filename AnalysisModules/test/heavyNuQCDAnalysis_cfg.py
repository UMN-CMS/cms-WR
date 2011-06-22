import FWCore.ParameterSet.Config as cms

import os

import PhysicsTools.PythonAnalysis.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

#import sys
#isMC=sys.modules['__main__'].isMC
#isMCsignal=sys.modules['__main__'].isMCsignal
#process = sys.modules['__main__'].process

isMC=False
isMCsignal=False
Training=False
isRun2010LoLumi=False
isRun2011=True

isData=not isMC

## Low and high lumi data selection is controlled by the JSON-derived cfi's imported
## below. For run 2010, the low lumi data is that for which the HLT_Mu9 trigger path
## was active and unprescaled, (uncertified) run range 133446 - 147116. Certification
## restricts this run range further.
##
if not isRun2011:
    isRun2010HiLumi=not isRun2010LoLumi
else:
    isRun2010HiLumi = False 
    isRun2010LoLumi = False 

process = cms.Process("PAT");

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
# source
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring( "file:foo.root" )
)
# process.load('HeavyNu.AnalysisModules.in_cff')

if isData:
    if isRun2011:
        print "===========> Flag is SET for 2011 luminosity data <============"
        from HeavyNu.AnalysisModules.goodLumiList_may10rereco_2011_mu24x_cfi import lumisToProcess
    else:
        if isRun2010LoLumi:
            print "===========> Flag is SET for 2010 LOW luminosity data <============"
            from HeavyNu.AnalysisModules.goodLumiList_apr21rereco_2010_mu9_cfi import lumisToProcess
        else:
            print "===========> Flag is SET for 2010 HIGH luminosity data <============"
            from HeavyNu.AnalysisModules.goodLumiList_apr21rereco_2010_mu15_cfi import lumisToProcess    
    process.source.lumisToProcess = lumisToProcess

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
if (isMC):
    print "=================> MC flag is SET <===================="
    process.GlobalTag.globaltag = cms.string('START41_V0::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V14::All')

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
# PAT Jet Energy Corrections - MC vs Data
########################################
# 

from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection
if isMC:
    switchJetCollection( process,
                         jetCollection=cms.InputTag('ak5CaloJets'),
                         jetCorrLabel=('AK5Calo', ['L1Offset','L2Relative','L3Absolute']))
else:
    switchJetCollection( process,
                         jetCollection=cms.InputTag('ak5CaloJets'),
                         jetCorrLabel=('AK5Calo', ['L1Offset','L2Relative','L3Absolute']))

########################################
# PAT Trigger matching
########################################
# imported directly from PhysicsTools/PatExamples/test/analyzePatTrigger_onTheFly_cfg.py
#
process.load("HeavyNu.AnalysisModules.hnutrigmatch_cfi")

### ============
### Python tools
### ============
### Attention: order matters!

## --
## Switch to selected PAT objects in the main work flow
## --
from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
removeCleaning( process, isData )

## Special change for saving good products in data
if isData:
    process.out.outputCommands = cms.untracked.vstring("keep *")

## --
## Switch on PAT trigger - but only for data!
## --
from PhysicsTools.PatAlgos.tools.trigTools import *
if isData:
    switchOnTriggerMatching( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ] )
    removeCleaningFromTriggerMatching( process )
    if isRun2011:
        process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu24_v*" )' )
    else:
        if isRun2010LoLumi: process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu9" )' )
        else:               process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu15_v1" )' )

##########################################
## Add analysis
##########################################

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("anal.root"),
)

process.hNuCalo = cms.EDFilter("MuJetBackground",
    trigMatchPset = cms.PSet(
        trigEventTag = cms.InputTag( "" ),
        muonMatch    = cms.string( '' ),
        randomSeed   = cms.int32( 0 ),  # for MC
        year         = cms.int32( 2011 )  # for MC
    ),
    muIDPset = cms.PSet(
        eraForId     = cms.int32( 2011 )
    ),
    pileupEra         = cms.int32(20110),
    DoLog        = cms.bool( False ),
    muonTag      = cms.InputTag( 'selectedPatMuons' ),
    jetTag       = cms.InputTag( 'selectedPatJets' ),
    electronTag  = cms.InputTag( 'selectedPatElectrons' ),
    photonTag    = cms.InputTag( 'selectedPatPhotons' ),
    BtagName     = cms.string('jetProbabilityBJetTags'),
    minMu1pt     = cms.double(60.),
    minMu2pt     = cms.double(30.),
    minJetPt     = cms.double(40),
    maxMuAbsEta  = cms.double(2.4),
    maxJetAbsEta = cms.double(2.5),
    minMuMuMass  = cms.double(200),
    min4objMass  = cms.double(520),
    minMuonJetdR = cms.double(0.5), 
    
    # muonTrackIsoLimitGeV  = cms.double(10.0),
    muonTrackIsoLimitGeV  = cms.double(0.10),
    dimuonMaxVertexZsepCM = cms.double(0.03),
    
    ZmassWinMinGeV= cms.double(84.),
    ZmassWinMaxGeV= cms.double(98.),

    applyJECUsign  = cms.int32(0),
    applyMESfactor = cms.double(1.0),

    isSignal     = cms.bool(False),
    mNuRnormalization = cms.double(1000.0), 

    #--- Specific variables added for INR QCD cross-check ---#
    hybridSCs   = cms.InputTag( 'correctedHybridSuperClusters' ),
    multi5x5SCs = cms.InputTag( 'correctedMulti5x5SuperClustersWithPreshower' ),
    minimumSuperClusterEt = cms.double(10.0), 

    getSurvivalRate = cms.bool(False),    
    doClosureTest   = cms.bool(True),    
    doQuadJetTest   = cms.bool(True),    
    #--- New results from 42x re-reco 2010 (36/pb) data ---#
    # reweightPtLow  = cms.vdouble( 20,25,30,40,60,100 ),
    # reweightPtHigh = cms.vdouble( 25,30,40,60,100,1000 ),
    # reweightLoose  = cms.vdouble( 0.0493001,0.0512686,0.0625999,0.0822622,0.159119,0.273381 ),
    # reweightTight  = cms.vdouble( 0.0609952,0.0597178,0.0608489,0.062275,0.0955056,0.157895 ),
    #--- New results from 42x re-reco 2011 (204/pb) data ---#
    reweightPtLow  = cms.vdouble( 30,40,60,100 ),
    reweightPtHigh = cms.vdouble( 40,60,100,1000 ),
    reweightLoose  = cms.vdouble( 0.0652302,0.0854414,0.147673,0.311429 ),
    reweightTight  = cms.vdouble( 0.0656517,0.0671176,0.0775541,0.177215 ),

    minimumMuJetdPhi = cms.double(2.8274334),
    minimumJetPtForDijets = cms.double(10.),
    minimumDeltaRforExtraJets = cms.double(0.7),
    METvariety = cms.int32(1) 
)

# process.load("HeavyNu.AnalysisModules.heavynuanalysis_cfi")

if isData:
    # turn on trigger match requirement
    process.hNuCalo.trigMatchPset.trigEventTag=cms.InputTag("patTriggerEvent")
    process.hNuCalo.trigMatchPset.muonMatch=cms.string('muonTriggerMatchHLTMuons')
else:
    # turn on MC trigger simulation
    process.hNuCalo.trigMatchPset.randomSeed=cms.int32(os.getpid())

if Training:
    process.hNuCalo.trainingFileName=cms.untracked.string("changeme_nntraining.txt")
    
if isMCsignal:
    process.hNuCalo.isSignal = cms.bool(True)

if isMC:
    process.hNuCalo2010 = process.hNuCalo.clone(minMu2pt = cms.double(20.))
    process.hNuCalo2010.pileupEra          = cms.int32(20100)
    process.hNuCalo2010.muIDPset.eraForId  = cms.int32(2010)
    process.hNuCalo2010.trigMatchPset.year = cms.int32(2010)
    process.hNuCalo2010.reweightPtLow  = cms.vdouble( 20,25,30,40,60,100 ),
    process.hNuCalo2010.reweightPtHigh = cms.vdouble( 25,30,40,60,100,1000 ),
    process.hNuCalo2010.reweightLoose  = cms.vdouble( 0.0493001,0.0512686,0.0625999,0.0822622,0.159119,0.273381 ),
    process.hNuCalo2010.reweightTight  = cms.vdouble( 0.0609952,0.0597178,0.0608489,0.062275,0.0955056,0.157895 ),

    process.p += process.hNuCalo2010
    process.p1 = cms.Path( process.patDefaultSequence + process.hNuCalo )
    process.s  = cms.Schedule(process.p,process.p1)
)

else:
    process.p += process.hNuCalo
