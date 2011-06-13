import FWCore.ParameterSet.Config as cms

import os

#import sys
#isMC=sys.modules['__main__'].isMC
#isMCsignal=sys.modules['__main__'].isMCsignal
#process = sys.modules['__main__'].process

isMC=True
isMCsignal=False
Training=False
isRun2010LoLumi=False
isRun2011=True
isPileupMC=True

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
    fileNames=cms.untracked.vstring('file:input.root')
)

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
    process.GlobalTag.globaltag = cms.string('GR_R_42_V14::All')

process.load("Configuration.StandardSequences.MagneticField_cff")

################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

# to activate new corrections
## process.load("RecoEgamma.EgammaElectronProducers.correctedGsfElectrons_cfi")
## process.correctedGsfElectrons.applyEtaCorrection = cms.bool(True)

## pat sequences to be loaded:
process.load("PhysicsTools.PatAlgos.patSequences_cff")
## process.makePatElectrons = cms.Sequence(process.correctedGsfElectrons * process.electronMatch * process.patElectrons) 

########################################
# Output module - has to be defined before PAT python tools will work
########################################

if isData:
    from PhysicsTools.PatAlgos.tools.coreTools import *
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

## --
## Switch on PAT trigger - but only for data!
## --
from PhysicsTools.PatAlgos.tools.trigTools import *
if isData:
    switchOnTriggerMatching( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ] )
    removeCleaningFromTriggerMatching( process, outputModule = '' )
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

process.load("HeavyNu.AnalysisModules.heavynutopanalysis_cfi")
process.hNuTop.minMu2pt         = cms.double(30.)
process.hNuTop.studyScaleFactor = cms.bool(True)

if isMC:
    process.hNuTop.applyMuIDEffcorr = cms.bool(True)
    process.hNuTop.applyEleEScale = cms.bool(False) 
    process.hNuTop.applyEleIDweight = cms.bool(True) 
    #--- 2010 values ---#
    # process.hNuTop.EBidWgt = cms.double( 0.978 ) 
    # process.hNuTop.EEidWgt = cms.double( 0.994 ) 
    #--- 2011 values ---#
    process.hNuTop.EBidWgt = cms.double( 0.975 ) 
    process.hNuTop.EEidWgt = cms.double( 0.964 ) 
else:
    process.hNuTop.applyMuIDEffcorr = cms.bool(False)
    process.hNuTop.applyEleIDweight = cms.bool(False) 
    process.hNuTop.applyEleEScale = cms.bool(False) 
    # process.hNuTop.EBscalefactor = cms.double( 1.0 ) 
    # process.hNuTop.EEscalefactor = cms.double( 1.04 ) 

if isData:
    # turn on trigger match requirement
    process.hNuTop.trigMatchPset.trigEventTag=cms.InputTag("patTriggerEvent")
    process.hNuTop.trigMatchPset.muonMatch=cms.string('muonTriggerMatchHLTMuons')
else:
    # turn on MC trigger simulation
    process.hNuTop.trigMatchPset.randomSeed=cms.int32(os.getpid())

if Training:
    process.hNuTop.trainingFileName=cms.untracked.string("changeme_nntraining.txt")

#--- Switch to correctedGsfElectrons ---#
#process.eidCutBasedExt.src = cms.InputTag("correctedGsfElectrons")
#process.eidRobustLoose.src = cms.InputTag("correctedGsfElectrons")
#process.eidRobustTight.src = cms.InputTag("correctedGsfElectrons")
#process.eidRobustHighEnergy.src = cms.InputTag("correctedGsfElectrons")
#process.eidLoose.src = cms.InputTag("correctedGsfElectrons")
#process.eidTight.src = cms.InputTag("correctedGsfElectrons")
#process.eleIsoDepositEcalFromHits.src = cms.InputTag("correctedGsfElectrons")
#process.eleIsoDepositHcalFromTowers.src = cms.InputTag("correctedGsfElectrons")
#process.eleIsoDepositTk.src = cms.InputTag("correctedGsfElectrons")
#process.electronMatch.src = cms.InputTag("correctedGsfElectrons")
#process.softElectronCands.electrons = cms.InputTag("correctedGsfElectrons")
#process.softElectronSelector.input = cms.InputTag("correctedGsfElectrons")
#process.softElectronTagInfos.leptons = cms.InputTag("correctedGsfElectrons")
#process.patElectrons.electronSource = cms.InputTag("correctedGsfElectrons")
    
## ---
## Define the paths
## ---
if isMC:
    process.hNu2011Top20 = process.hNuTop.clone(minMu2pt = cms.double(20.))
    process.hNu2011Top30 = process.hNuTop.clone(minMu2pt = cms.double(30.))
    process.hNu2010Top20 = process.hNuTop.clone(minMu2pt = cms.double(20.))
    process.hNu2010Top30 = process.hNuTop.clone(minMu2pt = cms.double(30.))
    process.hNu2010Top20.pileupEra = cms.int32(20100)
    process.hNu2010Top30.pileupEra = cms.int32(20100)
    process.hNu2010Top20.muIDPset.eraForId  = cms.int32(2010)
    process.hNu2010Top30.muIDPset.eraForId  = cms.int32(2010)
    process.hNu2010Top20.trigMatchPset.year = cms.int32(2010)
    process.hNu2010Top30.trigMatchPset.year = cms.int32(2010)
    process.hNu2010Top20.EBidWgt = cms.double( 0.978 ) 
    process.hNu2010Top30.EBidWgt = cms.double( 0.978 ) 
    process.hNu2010Top20.EEidWgt = cms.double( 0.994 ) 
    process.hNu2010Top30.EEidWgt = cms.double( 0.994 ) 
else:
    process.hNuTop20 = process.hNuTop.clone(minMu2pt = cms.double(20.))
    process.hNuTop30 = process.hNuTop.clone(minMu2pt = cms.double(30.))

if isMC:
    process.p2010Top20 = cms.Path(process.patDefaultSequence+process.hNu2010Top20)
    process.p2010Top30 = cms.Path(process.patDefaultSequence+process.hNu2010Top30)
    process.p2011Top20 = cms.Path(process.patDefaultSequence+process.hNu2011Top20)
    process.p2011Top30 = cms.Path(process.patDefaultSequence+process.hNu2011Top30)
    process.s = cms.Schedule(process.p2010Top20,
                             process.p2010Top30,
                             process.p2011Top20,
                             process.p2011Top30)
else:
    process.pTop20 = cms.Path(process.patDefaultSequence+process.hNuTop20)
    process.pTop30 = cms.Path(process.patDefaultSequence+process.hNuTop30)
    process.s = cms.Schedule(process.pTop20,
                             process.pTop30)
