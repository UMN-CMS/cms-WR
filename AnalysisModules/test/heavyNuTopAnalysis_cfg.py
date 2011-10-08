import FWCore.ParameterSet.Config as cms

import os

#import sys
#isMC=sys.modules['__main__'].isMC
#isMCsignal=sys.modules['__main__'].isMCsignal
#process = sys.modules['__main__'].process

isMC=False
isMCsignal=False
Training=False
isRun2011LoLumi=False
isRun2011VeryHiLumi=True
isPileupMC=True
# Note: Assumes running PF.  Do not change isPFJets to False!
isPFJets=True

isData=not isMC

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
    fileNames=cms.untracked.vstring('/store/data/Run2011A/SingleMu/AOD/05Aug2011-v1/0000/96291883-D9C0-E011-9352-0026189437F8.root')
)

if isData:
    if isRun2011LoLumi:
        print "===========> Flag is SET for 2011 LOW luminosity data <============"
        from HeavyNu.AnalysisModules.goodLumiList_160431_163869_Mu24_cfi import lumisToProcess
    else:
        if isRun2011VeryHiLumi:
            print "===========> Flag is SET for 2011 HIGH luminosity data <============"
            from HeavyNu.AnalysisModules.goodLumiList_173236_Mu40eta2p1_cfi import lumisToProcess
        else:
            print "===========> Flag is SET for 2011 MEDIUM luminosity data <============"
            from HeavyNu.AnalysisModules.goodLumiList_165088_173198_Mu40_cfi import lumisToProcess

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
        print "=============> isPileupMC flag is SET <================"
        process.GlobalTag.globaltag=cms.string('START42_V13::All')
else:
    print "===============> Running on DATA <===================="
    process.GlobalTag.globaltag = cms.string('GR_R_42_V20::All')

process.load("Configuration.StandardSequences.MagneticField_cff")

################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

# to activate new corrections
## process.load("RecoEgamma.EgammaElectronProducers.correctedGsfElectrons_cfi")
## process.correctedGsfElectrons.applyEtaCorrection = cms.bool(True)

## pat sequences to be loaded:
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("RecoMuon.MuonIdentification.refitMuons_cfi")
## process.makePatElectrons = cms.Sequence(process.correctedGsfElectrons * process.electronMatch * process.patElectrons) 

from PhysicsTools.PatAlgos.tools.pfTools import *
########################################
# Output module - has to be defined before PAT python tools will work but is NOT used
########################################
process.out = cms.OutputModule("PoolOutputModule",
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

if isData:
    from PhysicsTools.PatAlgos.tools.coreTools import *
    removeMCMatchingPF2PAT(process, '')
    removeMCMatching(process, ['All'], outputInProcess = False)

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
        if isRun2011VeryHiLumi:
            process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu40_eta2p1_v*",1,0 )' )
        else:
            process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu40_v*",1,0 )' )

##########################################
## Add analysis
##########################################

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("anal.root"),
)

process.load("HeavyNu.AnalysisModules.heavynutopanalysis_cfi")
process.hNuTop.minMu2pt         = cms.double(30.)
process.hNuTop.studyScaleFactor = cms.bool(False)
process.hNuTop.jetTag           = cms.InputTag( 'selectedPatJetsPFlow' )
if isData:
    process.hNuTop.muonTag = cms.InputTag( 'selectedPatMuonsTriggerMatch' )

if isRun2011LoLumi:
    process.hNuTop.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu24_v1','HLT_Mu24_v2' )
    process.hNuTop.trigMatchPset.triggerPt = cms.double( 24. )
else:
    process.hNuTop.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_v1','HLT_Mu40_v2','HLT_Mu40_v3','HLT_Mu40_v5','HLT_Mu40_eta2p1_v1' ) 
    process.hNuTop.trigMatchPset.triggerPt = cms.double( 40. )

if isMC:
    # Do not know any of these values yet, so assume no correction
    process.hNuTop.applyMuIDEffcorr   = cms.bool(False)
    process.hNuTop.applyEleEScale     = cms.bool(False) 
    process.hNuTop.applyEleIDweight   = cms.bool(False) 
    process.hNuTop.pileupEra          = cms.int32(20110)
    process.hNuTop.muIDPset.eraForId  = cms.int32(2011)
    process.hNuTop.trigMatchPset.year = cms.int32(2011)
    process.hNuTop.EBidWgt = cms.double( 1.000 ) 
    process.hNuTop.EEidWgt = cms.double( 1.000 ) 
else:
    process.hNuTop.muonTag          = cms.InputTag( 'selectedPatMuonsTriggerMatch' )
    process.hNuTop.applyMuIDEffcorr = cms.bool(False)
    process.hNuTop.applyEleIDweight = cms.bool(False) 
    process.hNuTop.applyEleEScale   = cms.bool(False) 

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
    process.hNuTopLoLumi = process.hNuTop.clone()
    process.hNuTopLoLumi.trigMatchPset.triggerPt = cms.double(24.)

    process.hNuTopHiLumi = process.hNuTop.clone()
    process.hNuTopHiLumi.trigMatchPset.triggerPt = cms.double(40.)

    process.pLo = cms.Path( process.patDefaultSequence * process.refitMuons * process.modifiedPF2PATSequence * process.hNuTopLoLumi ) 
    process.pHi = cms.Path( process.patDefaultSequence * process.refitMuons * process.modifiedPF2PATSequence * process.hNuTopHiLumi ) 
    process.s   = cms.Schedule( process.pLo,process.pHi )

else:
    process.p = cms.Path( process.patDefaultSequence * process.refitMuons * process.modifiedPF2PATSequence * process.hNuTop ) 

