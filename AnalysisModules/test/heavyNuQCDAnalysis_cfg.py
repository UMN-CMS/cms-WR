import FWCore.ParameterSet.Config as cms

import os

import PhysicsTools.PythonAnalysis.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

isMC=False
isMCsignal=False
Training=False
isRun2011LoLumi=True
isDijetStudy=True
isPFJets=True

isData=not isMC

## Low and high lumi data selection is controlled by the JSON-derived cfi's imported
## below. For run 2010, the low lumi data is that for which the HLT_Mu9 trigger path
## was active and unprescaled, (uncertified) run range 133446 - 147116. Certification
## restricts this run range further.
##
isRun2011HiLumi = not isRun2011LoLumi

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

if isData:
#    if isDijetStudy:
    from HeavyNu.AnalysisModules.goodLumiList_160404_176309_JSON_cfi import lumisToProcess
#    else:
#        if isRun2011LoLumi:
#            print "===========> Flag is SET for 2011 LOW luminosity data <============"
#            from HeavyNu.AnalysisModules.goodLumiList_160431_163869_Mu24_cfi.py import lumisToProcess
#        else:
#            print "===========> Flag is SET for 2011 HIGH luminosity data <============"
#            from HeavyNu.AnalysisModules. import lumisToProcess    
    process.source.lumisToProcess = lumisToProcess

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
if (isMC):
    print "=================> MC flag is SET <===================="
    process.GlobalTag.globaltag=cms.string('START42_V13::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V20::All')

process.load("Configuration.StandardSequences.MagneticField_cff")

################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

## pat sequences to be loaded:
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("RecoMuon.MuonIdentification.refitMuons_cfi")
from PhysicsTools.PatAlgos.tools.pfTools import *

########################################
# Output module - has to be defined before PAT python tools will work
########################################
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('heavynu_candevents.root'),
                               SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('p')),
                               outputCommands = cms.untracked.vstring("keep *")
                              )

if isPFJets:
    postfix = "PFlow"
    if isMC:
        usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=isMC, postfix=postfix, 
                  jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute']))
    else:
        usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=isMC, postfix=postfix, 
                  jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']))
    # Remove pileup, muon, and electron candidates from jets 
    # N.B.: This should already be done by default
    getattr(process,"pfNoPileUp"+postfix).enable   = True
    getattr(process,"pfNoMuon"+postfix).enable     = True
    getattr(process,"pfNoElectron"+postfix).enable = True
    # turn to false when running on data
    getattr(process, "patElectrons"+postfix).embedGenMatch = isMC
    getattr(process, "patMuons"+postfix).embedGenMatch     = isMC
    getattr(process, "patTaus"+postfix).embedGenMatch      = isMC
    # Needed for jet energy corrections
    process.pfPileUpPFlow.Enable = True
    process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)
    process.pfPileUpPFlow.Vertices = cms.InputTag('goodOfflinePrimaryVertices')
    process.pfJetsPFlow.doAreaFastjet = True
    process.pfJetsPFlow.doRhoFastjet = False

    #-----------------------------------------------------------------------#
    #--- Jet Corrections using PF and PF2PAT                             ---#
    #--- twiki reference: CMSPublic/WorkBookJetEnergyCorrections         ---#
    #--- See also: PhysicsTools/PatExamples/test/patTuple_42x_jec_cfg.py ---#
    #-----------------------------------------------------------------------#
    from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
    process.goodOfflinePrimaryVertices = cms.EDFilter(
        "PrimaryVertexObjectFilter",
        filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
        src=cms.InputTag('offlinePrimaryVertices')
    )

    # Compute the mean pt per unit area (rho) from the PFchs inputs
    from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
    process.kt6PFJetsPFlow = kt4PFJets.clone(
        rParam = cms.double(0.6),
        src = cms.InputTag('pfNoElectron'+postfix),
        doAreaFastjet = cms.bool(True),
        doRhoFastjet = cms.bool(True)
    )
    process.patJetCorrFactorsPFlow.rho = cms.InputTag("kt6PFJetsPFlow", "rho")

    # Add the PV selector and KT6 producer to the sequence
    getattr(process,"patPF2PATSequence"+postfix).replace(
        getattr(process,"pfNoElectron"+postfix),
        getattr(process,"pfNoElectron"+postfix)*process.kt6PFJetsPFlow 
    )

    process.modifiedPF2PATSequence = cms.Sequence(    
        process.goodOfflinePrimaryVertices*
        getattr(process,"patPF2PATSequence"+postfix)
    )

## ---
## Define the path
## ---
process.p = cms.Path(
  process.patDefaultSequence * process.refitMuons
)
if isPFJets:
    process.p += process.modifiedPF2PATSequence

if isData:
    # process.outpath  = cms.EndPath(process.out)
    from PhysicsTools.PatAlgos.tools.coreTools import *
    removeMCMatchingPF2PAT( process, '' )
    removeMCMatching(process, ['All'], outputInProcess = False)

########################################
# PAT Jet Energy Corrections - MC vs Data
########################################
# 

# from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection
# if isMC:
#     switchJetCollection( process,
#                          jetCollection=cms.InputTag('ak5CaloJets'),
#                          jetCorrLabel=('AK5Calo', ['L1Offset','L2Relative','L3Absolute']))
# else:
#     switchJetCollection( process,
#                          jetCollection=cms.InputTag('ak5CaloJets'),
#                          jetCorrLabel=('AK5Calo', ['L1Offset','L2Relative','L3Absolute']))
process.patJetCorrFactors.useRho = cms.bool(False)

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
removeCleaning( process, False )

## Special change for saving good products in data
if isData:
    process.out.outputCommands = cms.untracked.vstring("keep *")

## --
## Switch on PAT trigger - but only for data!
## --
from PhysicsTools.PatAlgos.tools.trigTools import *
if isData:
    switchOnTriggerMatching( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ], outputModule = '' )
    switchOnTriggerMatchEmbedding( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ], outputModule = '' )
    removeCleaningFromTriggerMatching( process, outputModule = '' )
    if isRun2011LoLumi:
        process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu24_v*",1,0 )' )
    else:
        process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu40_*",1,0 )' )

##########################################
## Add analysis
##########################################

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("anal.root"),
)

process.hNuQCD = cms.EDFilter("MuJetBackground",
    trigMatchPset = cms.PSet(
        trigEventTag = cms.InputTag( "" ),
        muonMatch    = cms.string( 'muonTriggerMatchHLTMuons' ),
        muonTriggers = cms.vstring( 'HLT_Mu24_v1','HLT_Mu24_v2','HLT_Mu40_v1','HLT_Mu40_v2','HLT_Mu40_v3','HLT_Mu40_v5','HLT_Mu40_eta2p1_v1' ),
        triggerPt    = cms.double( 40. ),
        firstRun     = cms.vint32( 0 ),
        lastRun      = cms.vint32( 999999 ),
        randomSeed   = cms.int32( 0 ),  # for MC
        year         = cms.int32( 2011 )  # for MC
    ),
    muIDPset = cms.PSet(
        eraForId     = cms.int32( 2011 )
    ),
    pileupEra         = cms.int32(20110),
    DoLog        = cms.bool( False ),
    isPFJets     = cms.bool( True ), 
    muonTag      = cms.InputTag( 'selectedPatMuons' ),
    jetTag       = cms.InputTag( 'selectedPatJetsPFlow' ),
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
    maxJetVZsepCM         = cms.double(0.1),

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

    getSurvivalRate = cms.bool(True),    
    doClosureTest   = cms.bool(False),    
    doQuadJetTest   = cms.bool(False),    
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
    # Take PF MET 
    METvariety = cms.int32(2) 
)

# process.load("HeavyNu.AnalysisModules.heavynuanalysis_cfi")

if isData:
    # turn on trigger match requirement
    process.hNuQCD.trigMatchPset.trigEventTag=cms.InputTag("patTriggerEvent")
    process.hNuQCD.trigMatchPset.muonMatch=cms.string('muonTriggerMatchHLTMuons')
else:
    # turn on MC trigger simulation
    process.hNuQCD.trigMatchPset.randomSeed=cms.int32(os.getpid())

if Training:
    process.hNuQCD.trainingFileName=cms.untracked.string("changeme_nntraining.txt")
    
if isMCsignal:
    process.hNuQCD.isSignal = cms.bool(True)

if isMC:
    process.hNuQCD2010 = process.hNuQCD.clone(minMu2pt = cms.double(20.))
    process.hNuQCD2010.pileupEra          = cms.int32(20100)
    process.hNuQCD2010.muIDPset.eraForId  = cms.int32(2010)
    process.hNuQCD2010.trigMatchPset.year = cms.int32(2010)
    process.hNuQCD2010.reweightPtLow  = cms.vdouble( 20,25,30,40,60,100 ),
    process.hNuQCD2010.reweightPtHigh = cms.vdouble( 25,30,40,60,100,1000 ),
    process.hNuQCD2010.reweightLoose  = cms.vdouble( 0.0493001,0.0512686,0.0625999,0.0822622,0.159119,0.273381 ),
    process.hNuQCD2010.reweightTight  = cms.vdouble( 0.0609952,0.0597178,0.0608489,0.062275,0.0955056,0.157895 ),

    process.p += process.hNuQCD2010
    process.p1 = cms.Path( process.patDefaultSequence + process.hNuQCD )
    process.s  = cms.Schedule(process.p,process.p1)

else:
    process.p += process.hNuQCD
