import FWCore.ParameterSet.Config as cms

import os

#import sys
#isMC=sys.modules['__main__'].isMC
#isMCsignal=sys.modules['__main__'].isMCsignal
#process = sys.modules['__main__'].process

isMC=False
isMCsignal=False
Training=False
isRun2010LoLumi=True
isRun2011=False

isData=not isMC

## Low and high lumi data selection is controlled by the JSON-derived cfi's imported
## below. For run 2010, the low lumi data is that for which the HLT_Mu9 trigger path
## was active and unprescaled, (uncertified) run range 133446 - 147116. Certification
## restricts this run range further.
##
isRun2010HiLumi=not isRun2010LoLumi

process = cms.Process("PAT");

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
# source
#process.source = cms.Source("PoolSource",
#    fileNames=cms.untracked.vstring('file:/hdfs/cms/skim/mu/39X/Dec22ReReco/Run2010B/HiLumi/pool_1_1_4Cv.root')
#    fileNames=cms.untracked.vstring('file:/local/cms/phedex/store/mc/Fall10/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/START38_V12-v3/0001/02D6EA60-9D02-E011-A091-90E6BA442F11.root')
#)
process.load('HeavyNu.AnalysisModules.in_cff')

if isData:
    if isRun2010LoLumi:
        print "===========> Flag is SET for LOW luminosity data <============"
    else:
        print "===========> Flag is SET for HIGH luminosity data <============"
    
    from HeavyNu.AnalysisModules.goodRunList_cfi import lumisToProcess
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

########################################
# Output module - has to be defined before PAT python tools will work
########################################

if isData:
<<<<<<< heavyNuTopAnalysis_cfg.py
#     process.out = cms.OutputModule("PoolOutputModule",
#                                    fileName = cms.untracked.string('heavynu_candevents.root'),
#     # save only events passing the full path
#                                    SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('p')),
#                                    outputCommands = cms.untracked.vstring("keep *")
#                                    )
#     process.outpath  = cms.EndPath(process.out)
      from PhysicsTools.PatAlgos.tools.coreTools import *
      removeMCMatching(process, ['All'], outputInProcess = False)
=======
    from PhysicsTools.PatAlgos.tools.coreTools import *
    removeMCMatching(process, ['All'], outputInProcess = False)
>>>>>>> 1.3

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
<<<<<<< heavyNuTopAnalysis_cfg.py
removeCleaning( process, False )

## Special change for saving good products in data
#if isData:
#    process.out.outputCommands = cms.untracked.vstring("keep *")
=======
removeCleaning( process, False )
>>>>>>> 1.3

## --
## Switch on PAT trigger - but only for data!
## --
from PhysicsTools.PatAlgos.tools.trigTools import *
if isData:
<<<<<<< heavyNuTopAnalysis_cfg.py
    switchOnTriggerMatching( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ], outputModule='' )
    removeCleaningFromTriggerMatching( process, outputModule='' )
=======
    switchOnTriggerMatching( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ], outputModule = '' )
    removeCleaningFromTriggerMatching( process, outputModule = '' )
>>>>>>> 1.3
    if isRun2010LoLumi: process.muonTriggerMatchHLTMuons.pathNames = cms.vstring('HLT_Mu9')
    else:               process.muonTriggerMatchHLTMuons.pathNames = cms.vstring('HLT_Mu15_v1')

##########################################
## Add analysis
##########################################

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("anal.root"),
)

process.load("HeavyNu.AnalysisModules.heavynutopanalysis_cfi")
process.hNuTop.minMu2pt = cms.double(30.)
if isMC:
    process.hNuTop.applyMuIDEffcorr = cms.bool(True)
    process.hNuTop.applyEleEScale = cms.bool(False) 
    process.hNuTop.applyEleIDweight = cms.bool(True) 
    process.hNuTop.EBidWgt = cms.double( 0.978 ) 
    process.hNuTop.EEidWgt = cms.double( 0.994 ) 
else:
    process.hNuTop.applyMuIDEffcorr = cms.bool(False)
    process.hNuTop.applyEleIDweight = cms.bool(False) 
    process.hNuTop.applyEleEScale = cms.bool(True) 
    process.hNuTop.EBscalefactor = cms.double( 1.0 ) 
    process.hNuTop.EEscalefactor = cms.double( 1.04 ) 

if isData:
    # turn on trigger match requirement
    process.hNuTop.trigMatchPset.trigEventTag=cms.InputTag("patTriggerEvent")
    process.hNuTop.trigMatchPset.muonMatch=cms.string('muonTriggerMatchHLTMuons')
else:
    # turn on MC trigger simulation
    process.hNuTop.trigMatchPset.randomSeed=cms.int32(os.getpid())

if Training:
    process.hNuTop.trainingFileName=cms.untracked.string("changeme_nntraining.txt")
    
## ---
## Define the paths
## ---
if isMC:
    process.hNu2011Top20 = process.hNuTop.clone(minMu2pt = cms.double(20.))
    process.hNu2011Top30 = process.hNuTop.clone(minMu2pt = cms.double(30.))
    process.hNu2010Top20 = process.hNuTop.clone(minMu2pt = cms.double(20.))
    process.hNu2010Top30 = process.hNuTop.clone(minMu2pt = cms.double(30.))
    process.hNu2010Top20.applyMuIDEffcorr   = cms.bool(False)
    process.hNu2010Top20.trigMatchPset.year = cms.int32(2010)
    process.hNu2010Top30.applyMuIDEffcorr   = cms.bool(False)
    process.hNu2010Top30.trigMatchPset.year = cms.int32(2010)
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
