import FWCore.ParameterSet.Config as cms

import os

#import sys
#isMC=sys.modules['__main__'].isMC
#isMCsignal=sys.modules['__main__'].isMCsignal
#process = sys.modules['__main__'].process

isMC=True
isMCsignal=False
Training=False
isRun2010LoLumi=True
isRun2011=False
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

## pat sequences to be loaded:
process.load("PhysicsTools.PatAlgos.patSequences_cff")

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
    switchOnTriggerMatching( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ], outputModule = '' )
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

process.load("HeavyNu.AnalysisModules.heavynuanalysis_cfi")
process.hNu.minMu2pt         = cms.double(30.)
process.hNu.studyScaleFactor = cms.bool(False)
process.hNu.studyMuSelectEff = cms.bool(True)

if isMC:
    process.hNu.applyMuIDEffcorr = cms.bool(True) # applied for both years now
    process.hNu.trigMatchPset.randomSeed=cms.int32(os.getpid())

if Training:
    process.hNu.trainingFileName=cms.untracked.string("changeme_nntraining.txt")
    
## ---
## Define the paths
## ---
if isMC:
    process.hNu2010 = process.hNu.clone(minMu2pt = cms.double(20.))
    process.hNu2010.pileupEra          = cms.int32(20100)
    process.hNu2010.muIDPset.eraForId  = cms.int32(2010)
    process.hNu2010.trigMatchPset.year = cms.int32(2010)
    process.hNu2011 = process.hNu.clone(minMu2pt = cms.double(30.))
    process.hNu2011.pileupEra          = cms.int32(20110)
    process.hNu2011.muIDPset.eraForId  = cms.int32(2011)
    process.hNu2011.trigMatchPset.year = cms.int32(2011)

    process.hNu2010MESHi = process.hNu2010.clone(applyMESfactor = cms.double(1.01))
    process.hNu2010MESLo = process.hNu2010.clone(applyMESfactor = cms.double(0.99))
    process.hNu2011MESHi = process.hNu2011.clone(applyMESfactor = cms.double(1.01))
    process.hNu2011MESLo = process.hNu2011.clone(applyMESfactor = cms.double(0.99))

    process.hNu2010TrigHi = process.hNu2010.clone(applyTrigEffsign = cms.int32(1))
    process.hNu2010TrigLo = process.hNu2010.clone(applyTrigEffsign = cms.int32(-1))
    process.hNu2011TrigHi = process.hNu2011.clone(applyTrigEffsign = cms.int32(1))
    process.hNu2011TrigLo = process.hNu2011.clone(applyTrigEffsign = cms.int32(-1))

    process.hNu2010IDHi = process.hNu2010.clone(applyMuIDEffsign = cms.int32(1))
    process.hNu2010IDLo = process.hNu2010.clone(applyMuIDEffsign = cms.int32(-1))
    process.hNu2011IDHi = process.hNu2011.clone(applyMuIDEffsign = cms.int32(1))
    process.hNu2011IDLo = process.hNu2011.clone(applyMuIDEffsign = cms.int32(-1))

if isMC:
    process.p2010       = cms.Path(process.patDefaultSequence+process.hNu2010)
    process.p2010MESHi  = cms.Path(process.patDefaultSequence+process.hNu2010MESHi)
    process.p2010MESLo  = cms.Path(process.patDefaultSequence+process.hNu2010MESLo)
    process.p2010TrigHi = cms.Path(process.patDefaultSequence+process.hNu2010TrigHi)
    process.p2010TrigLo = cms.Path(process.patDefaultSequence+process.hNu2010TrigLo)
    process.p2010IDHi   = cms.Path(process.patDefaultSequence+process.hNu2010IDHi)
    process.p2010IDLo   = cms.Path(process.patDefaultSequence+process.hNu2010IDLo)

    process.p2011       = cms.Path(process.patDefaultSequence+process.hNu2011)
    process.p2011MESHi  = cms.Path(process.patDefaultSequence+process.hNu2011MESHi)
    process.p2011MESLo  = cms.Path(process.patDefaultSequence+process.hNu2011MESLo)
    process.p2011TrigHi = cms.Path(process.patDefaultSequence+process.hNu2011TrigHi)
    process.p2011TrigLo = cms.Path(process.patDefaultSequence+process.hNu2011TrigLo)
    process.p2011IDHi   = cms.Path(process.patDefaultSequence+process.hNu2011IDHi)
    process.p2011IDLo   = cms.Path(process.patDefaultSequence+process.hNu2011IDLo)

    process.s = cms.Schedule(process.p2010,
                             process.p2010MESHi,
                             process.p2010MESLo,
                             process.p2010TrigHi,
                             process.p2010TrigLo,
                             process.p2010IDHi,
                             process.p2010IDLo,
                             process.p2011,
                             process.p2011MESHi,
                             process.p2011MESLo,
                             process.p2011TrigHi,
                             process.p2011TrigLo,
                             process.p2011IDHi,
                             process.p2011IDLo)
