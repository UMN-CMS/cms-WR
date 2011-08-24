import FWCore.ParameterSet.Config as cms

import os

#import sys
#isMC=sys.modules['__main__'].isMC
#isMCsignal=sys.modules['__main__'].isMCsignal
#process = sys.modules['__main__'].process

isMC=False
isMCsignal=False
Training=False
isRun2011LoLumi=True
# isRun2011=True
isPileupMC=True
isPFJets=False

isData=not isMC

#--- Things to change ---#
#--- Check the trigger match ---#
#--- Check the Global Tag ---#
#--- Check the JSON ---#
#--- Check the name of the module ---#

## Low and high lumi data selection is controlled by the JSON-derived cfi's imported
## below. For run 2010, the low lumi data is that for which the HLT_Mu9 trigger path
## was active and unprescaled, (uncertified) run range 133446 - 147116. Certification
## restricts this run range further.
##
isRun2011HiLumi=not isRun2011LoLumi

process = cms.Process("PAT");

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)

# source
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring('file:/local/cms/user/jmmans/heavyNuS11/HeavyNu_S11_AODSIM_1000_600/HeavyNu_S11_AODSIM_1000_600_001.root')
)

if isData:
    if isRun2011LoLumi or isRun2011HiLumi:
        print "===========> Flag is SET for 2011 luminosity data <============"
        # from HeavyNu.AnalysisModules.goodLumiList_mu24v1_cfi import lumisToProcess
        # from HeavyNu.AnalysisModules.goodLumiList_mu24v2_cfi import lumisToProcess
        # from HeavyNu.AnalysisModules.goodLumiList_mu40v1v2_cfi import lumisToProcess
        # from HeavyNu.AnalysisModules.goodLumiList_mu40v3_cfi import lumisToProcess
        from HeavyNu.AnalysisModules.goodLumiList_mu40v5_cfi import lumisToProcess
    # else:
    #     if isRun2010LoLumi:
    #         print "===========> Flag is SET for 2010 LOW luminosity data <============"
    #         from HeavyNu.AnalysisModules.goodLumiList_apr21rereco_2010_mu9_cfi import lumisToProcess
    #     else:
    #         print "===========> Flag is SET for 2010 HIGH luminosity data <============"
    #         from HeavyNu.AnalysisModules.goodLumiList_apr21rereco_2010_mu15_cfi import lumisToProcess    

    process.source.lumisToProcess = lumisToProcess

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
#
# Use 38X for MC, 39X for data
#
if (isMC):
    print "=================> MC flag is SET <===================="
    if (isPileupMC):
        process.GlobalTag.globaltag=cms.string('START41_V0::All')
        print "=============> isPileupMC flag is SET <================"
    else:
        print "========> Fall10 MC with Spring10 JEC applied <========"
        process.GlobalTag.globaltag = cms.string('START38_V14::All')
else:
    print "===============> Running on DATA <===================="
    # process.GlobalTag.globaltag = cms.string('GR_R_42_V20::All')
    process.GlobalTag.globaltag = cms.string('GR_P_V22::All')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')

################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

## pat sequences to be loaded:
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("RecoMuon.MuonIdentification.refitMuons_cfi")
# process.muonMatch.muonSource = cms.InputTag("refitMuons")
# process.patMuons.muonSource = cms.InputTag("refitMuons")
# process.makePatMuons = cms.Sequence( process.refitMuons * process.muonMatch * process.patMuons )

from PhysicsTools.PatAlgos.tools.pfTools import *

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('heavynu_candevents.root'),
                               # save only events passing the full path
                               SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('p')),
                               outputCommands = cms.untracked.vstring("keep *")
                               )
if isPFJets:
    postfix = "PFlow"
    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=isMC, postfix=postfix)
    # turn to false when running on data
    getattr(process, "patElectrons"+postfix).embedGenMatch = isMC
    getattr(process, "patMuons"+postfix).embedGenMatch = isMC
    getattr(process, "patTaus"+postfix).embedGenMatch = isMC
    
#------------------------------#
#--- Include Generic Tracks ---#
#------------------------------#
#--- Generic PAT tracks modules stolen from ElectroWeakAnalysis/Skimming/python ---#
process.patAODTrackCandsUnfiltered = cms.EDProducer("ConcreteChargedCandidateProducer",
    src          = cms.InputTag("generalTracks"),
    particleType = cms.string('mu+')   # to fix mass hypothesis
)
process.patAODTrackCands = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("patAODTrackCandsUnfiltered"),
    cut = cms.string('pt > 10')
)
from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import patGenericParticles
process.allPatTracks = patGenericParticles.clone(
    src = cms.InputTag("patAODTrackCands")
)
from PhysicsTools.PatAlgos.selectionLayer1.trackSelector_cfi import *
process.patTracksPt10 = selectedPatTracks.clone(
    cut = 'pt > 10.'
)
process.patTrackSequence = cms.Sequence( 
        process.patAODTrackCandsUnfiltered *
        process.patAODTrackCands *
        process.allPatTracks *
        process.patTracksPt10
)

## ---
## Define the path
## ---
process.p = cms.Path(
  process.patDefaultSequence * process.patTrackSequence * process.refitMuons
)

if isPFJets:
    process.p += getattr(process,"patPF2PATSequence"+postfix)

########################################
# Output module - has to be defined before PAT python tools will work
########################################

if isData:
    # process.outpath  = cms.EndPath(process.out)
    from PhysicsTools.PatAlgos.tools.coreTools import *
    if isPFJets:
        removeMCMatchingPF2PAT( process, '' )
    removeMCMatching(process, ['All'], outputInProcess = False)
        

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
# removeCleaning( process, isData )

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
        process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu40_v*",1,0 )' )
    # else:
        # if isRun2010LoLumi: process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu9" )' )
        # else:               process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu15_v1" )' )

##########################################
## Add analysis
##########################################

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("anal.root"),
)

process.load("HeavyNu.AnalysisModules.heavynuanalysis_cfi")
if isMCsignal:
	process.load("HeavyNu.AnalysisModules.heavyNuGenFilter_cfi")
	process.hNuGenFilter.keepIds = cms.vint32(2,)

process.hNu.studyMuSelectEff = cms.bool(False)
process.hNu.studyScaleFactor = cms.bool(False)

process.hNu.minMu2pt         = cms.double(30.)
process.hNu.pileupEra        = cms.int32(20110)
process.hNu.applyMuIDEffcorr = cms.bool(isMC)
process.hNu.muonTag          = cms.InputTag( 'selectedPatMuonsTriggerMatch' )

process.hNu.isPFJets = cms.bool(isPFJets)
if isPFJets:
    process.hNu.jetTag  = cms.InputTag( 'selectedPatJetsPFlow')
    process.hNu.muonTag = cms.InputTag( 'selectedPatMuonsTriggerMatch')
    
if isData:
    # turn on trigger match requirement
    process.hNu.trigMatchPset.trigEventTag=cms.InputTag("patTriggerEvent")
    process.hNu.trigMatchPset.muonMatch=cms.string('muonTriggerMatchHLTMuons')
else:
    # turn on MC trigger simulation
    process.hNu.trigMatchPset.randomSeed=cms.int32(os.getpid())

if Training:
    process.hNu.trainingFileName=cms.untracked.string("changeme_nntraining.txt")

#--- Changes necessary to sort out trigger complications ---#
process.hNuMu24v1   = process.hNu.clone() 
process.hNuMu24v2   = process.hNu.clone() 
process.hNuMu40v1v2 = process.hNu.clone() 
process.hNuMu40v3   = process.hNu.clone() 
process.hNuMu40v5   = process.hNu.clone() 

process.hNuMu24v1.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu24_v1' )
process.hNuMu24v2.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu24_v2' )
process.hNuMu40v1v2.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_v1','HLT_Mu40_v2' )
process.hNuMu40v3.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_v3' )
process.hNuMu40v5.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_v5' )

process.hNuMu24v1.trigMatchPset.triggerPt = cms.double( 24. )
process.hNuMu24v2.trigMatchPset.triggerPt = cms.double( 24. )
process.hNuMu40v1v2.trigMatchPset.triggerPt = cms.double( 40. )
process.hNuMu40v3.trigMatchPset.triggerPt = cms.double( 40. )
process.hNuMu40v5.trigMatchPset.triggerPt = cms.double( 40. )
    
if isMCsignal:
    process.hNu.isSignal = cms.bool(True)
    process.p += process.hNuGenFilter*process.hNu
else:
    if isMC: 
        process.p += process.hNu
    else:
        # process.p += process.hNuMu24v1
        # process.p += process.hNuMu24v2
        # process.p += process.hNuMu40v1v2
        # process.p += process.hNuMu40v3
        process.p += process.hNuMu40v5
