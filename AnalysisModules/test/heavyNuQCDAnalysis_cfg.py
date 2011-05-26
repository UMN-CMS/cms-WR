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
filterByJSON=True
filterByGoodRuns=True

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

if isData and filterByJSON:
    if isRun2011:
        print "===========> Flag is SET for 2011 luminosity data <============"
        if filterByGoodRuns:
            from HeavyNu.AnalysisModules.run2011LumiJSONapr29_cfi import lumisToProcess
            process.source.lumisToProcess = lumisToProcess
        else:
            process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
            # myDCSLumis = LumiList.LumiList(filename = '/home/phys/dahmes/Work/Physics/WR/CMSSW_4_1_3/src/HeavyNu/AnalysisModules/test/dcsOnly_160404-163387.json').getCMSSWString().split(',')
            myDCSLumis = LumiList.LumiList(filename = '/home/phys/dahmes/Work/Physics/WR/CMSSW_4_1_3/src/HeavyNu/AnalysisModules/test/dcsOnly_160404-163796.json').getCMSSWString().split(',')
            process.source.lumisToProcess.extend(myDCSLumis)
    else:
        if isRun2010LoLumi:
            print "===========> Flag is SET for 2010 LOW luminosity data <============"
            from HeavyNu.AnalysisModules.run2010loLumiRunList_cfi import lumisToProcess
        else:
            print "===========> Flag is SET for 2010 HIGH luminosity data <============"
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
    # process.GlobalTag.globaltag = cms.string('GR_R_38X_V15::All')
    process.GlobalTag.globaltag = cms.string('GR_P_V17::All')

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
                         jetCorrLabel=('AK5Calo', ['L2Relative','L3Absolute']))
else:
    switchJetCollection( process,
                         jetCollection=cms.InputTag('ak5CaloJets'),
                         jetCorrLabel=('AK5Calo', ['L2Relative','L3Absolute','L2L3Residual']))

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
    process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu24_v*" )' )
    process.muonTriggerMatchHLTMuons.maxDeltaR = cms.double( 0.2 )
    process.muonTriggerMatchHLTMuons.maxDPtRel = cms.double( 1.0 )
    # if isRun2010LoLumi: process.muonTriggerMatchHLTMuons.pathNames = cms.vstring('HLT_Mu9')
    # else:               process.muonTriggerMatchHLTMuons.pathNames = cms.vstring('HLT_Mu15_v1')

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
    #--- Old results from 2010 data ---#
    # reweightPtLow   = cms.vdouble( 0,10,15,20,25,30,40,60,100,200,400),
    # reweightPtHigh  = cms.vdouble( 10,15,20,25,30,40,60,100,200,400,1000),
    # reweightLoose   = cms.vdouble( 0.226356,0.104505,0.058554,0.0505105,0.0554867,0.0587291,0.0769231,0.152395,0.263158,0.5,1 ),
    # reweightTight   = cms.vdouble( 0.0581472,0.0357932,0.063638,0.0628168,0.0629979,0.0587147,0.0646372,0.0949227,0.135135,1,1 ),
    #--- New results, 43.4/pb from 2011 data ---# 
    reweightPtLow  = cms.vdouble( 20,25,30,40,60,100 ),
    reweightPtHigh = cms.vdouble( 25,30,40,60,100,1000 ),
    reweightLoose  = cms.vdouble( 0.0587241,0.0549809,0.0633958,0.0872151,0.165605,0.464286 ),
    reweightTight  = cms.vdouble( 0.0657958,0.0632302,0.0589719,0.0699013,0.07109,0.266667 ),

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

process.p += process.hNuCalo
