import FWCore.ParameterSet.Config as cms

from operator import isSequenceType
import os

#--- Analysis Mode ---#
#options are HNUMU, HNUE, TOP, QCD, CLO
#analysisMode = 'HNUE'

#--- Data/MC switch ---#
isMC=True
isData=not isMC

smode = 0 #0 for Z+Jets 1 for other

#--- Special flag for 44x/Fall11 ---#
is44x=False
is42x=False

#--- Signal MC flags ---#
isMCsignal=False
Training=False

#--- Trigger-based luminosity flags ---#
#--- only one should be True        ---#
#--- LoLumi     --> HLT_Mu24        ---#
#--- HiLumi     --> HLT_Mu40        ---#
#--- VeryHiLumi --> HLT_Mu40_eta2p1 ---#
isRun2011LoLumi     = False
isRun2011HiLumi     = False
isRun2011VeryHiLumi = False
isRun2012           = True

#--- Flags for data taking era ---#
isRun2011A = False

#--- Should always be True ---#
isPileupMC = True
isPFJets   = True

#--- HEEP ID for electrons ---#
heepVersion = -1

process = cms.Process("PATSKIM");

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)

# source
process.source = cms.Source("PoolSource",
   fileNames=cms.untracked.vstring("/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v1/0000/D602DE83-F18B-E111-A721-001A64789CF0.root")
)

if isData:
    if isRun2012:
        print "===========> Flag is SET for 2012 data <============"
        from HeavyNu.AnalysisModules.goodLumiList_2012_dynamic_cfi import lumisToProcess
    else:
        if isRun2011Mu24:
            print "===========> Flag is SET for 2011 LOW luminosity data <============"
            if is44x:
                from HeavyNu.AnalysisModules.goodLumiList_novRereco_160404_163869_mu24_cfi import lumisToProcess
            else:
                from HeavyNu.AnalysisModules.goodLumiList_160404_163869_may10rereco_Mu24_cfi import lumisToProcess
        else:
            if isRun2011Mu40eta2p1:
                print "===========> Flag is SET for 2011 HIGH luminosity data <============"
                if is44x:
                    from HeavyNu.AnalysisModules.goodLumiList_novRereco_173236_180252_mu40_eta2p1_cfi import lumisToProcess
                else:
                    from HeavyNu.AnalysisModules.goodLumiList_173236_180252_Mu40eta2p1_cfi import lumisToProcess
            elif isRun2011Mu40:
                print "===========> Flag is SET for 2011 MEDIUM luminosity data <============"
                if is44x:
                    from HeavyNu.AnalysisModules.goodLumiList_novRereco_165088_173198_mu40_cfi import lumisToProcess
                else:
                    from HeavyNu.AnalysisModules.goodLumiList_165088_173198_Mu40_cfi import lumisToProcess

    process.source.lumisToProcess = lumisToProcess

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

## Global Tags:
if (isMC):
    print "=================> MC flag is SET <===================="
    if (isPileupMC):
        if (is44x):
            process.GlobalTag.globaltag=cms.string('START44_V12::All')
        elif is42x:
            process.GlobalTag.globaltag=cms.string('START42_V13::All')
        else:
            process.GlobalTag.globaltag=cms.string('START52_V9::All')
        print "=============> isPileupMC flag is SET <================"
    # else:
    #     print "========> Fall10 MC with Spring10 JEC applied <========"
    #     process.GlobalTag.globaltag = cms.string('START38_V14::All')
else:
    print "===============> Running on DATA <===================="
    if (isRun2012):
        process.GlobalTag.globaltag = cms.string('GR_R_50_V13::All')
    else:
        if (is44x):
            process.GlobalTag.globaltag = cms.string('GR_R_44_V13::All')
        else:
            process.GlobalTag.globaltag = cms.string('GR_R_42_V20::All')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')

################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

## pat sequences to be loaded:
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("RecoMuon.MuonIdentification.refitMuons_cfi")
process.refitMuons.src = cms.InputTag("muons")
process.myRefitMuonSequence = cms.Sequence( process.refitMuons )
from PhysicsTools.PatAlgos.tools.pfTools import *

#--- Output module:
#--- Must be defined before PAT python tools will work
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('heavynu_candevents.root'),
                               # save only events passing the full path
                               SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('p')),
                               outputCommands = cms.untracked.vstring("keep *")
                               )
if isPFJets:
    postfix = "PFlow"
    if isMC:
        usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=isMC, postfix=postfix,
                  jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute']))
    else:
        if isRun2012:
            usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=isMC, postfix=postfix,
                      jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']))
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

## --------------------- ##
## Define the basic path ##
## --------------------- ##
if isMC:
   # Gen Level Energy balance filter to fix Pythia6 lhe interface bug
   process.load("HeavyNu.AnalysisModules.hnuTotalKinematicsFilter_cfi")
   process.AnalysisIntroSequence = cms.Sequence(
       process.hnuTotalKinematicsFilter * process.patDefaultSequence * process.patTrackSequence * process.myRefitMuonSequence
   )
else:
   process.AnalysisIntroSequence = cms.Sequence(
       process.patDefaultSequence * process.patTrackSequence * process.myRefitMuonSequence
   )
if isPFJets:
    process.AnalysisIntroSequence += process.modifiedPF2PATSequence

process.p = cms.Path(
    process.AnalysisIntroSequence
)


if isData:
    # process.outpath  = cms.EndPath(process.out)
    from PhysicsTools.PatAlgos.tools.coreTools import *
    if isPFJets:
        removeMCMatchingPF2PAT( process, '' )
    if (is42x):
        removeMCMatching(process, ['All'], outputInProcess = False)
    else:
        removeMCMatching(process, ['All'], outputModules = [])


#--- Calo Jet Energy Corrections: No longer used ---#
process.patJetCorrFactors.useRho = cms.bool(False)
# Corrections for Calo Jets: no longer used
# if isMC:
#     switchJetCollection( process,
#                          jetCollection=cms.InputTag('ak5CaloJets'),
#                          jetCorrLabel=('AK5Calo', ['L1Offset','L2Relative','L3Absolute']))
# else:
#     switchJetCollection( process,
#                          jetCollection=cms.InputTag('ak5CaloJets'),
#                          jetCorrLabel=('AK5Calo', ['L1Offset','L2Relative','L3Absolute','L2L3Residual']))



## ============================== ##
## Python tools --> Order matters ##
## ============================== ##

from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
if (is44x):
    removeCleaning( process, outputModules = [] )
else:
    removeCleaning( process, False )

# Special change for saving good products in data
if isData:
    process.out.outputCommands = cms.untracked.vstring("keep *")

#--- Trigger matching ---#
process.load("HeavyNu.AnalysisModules.hnutrigmatch_cfi")
from PhysicsTools.PatAlgos.tools.trigTools import *
if isData:
    if runMuonAnalysis:
        switchOnTriggerMatching( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ], outputModule = '' )
        switchOnTriggerMatchEmbedding( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ], outputModule = '' )
        removeCleaningFromTriggerMatching( process, outputModule = '' )
        if isRun2011Mu24:
            process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu24_v*",1,0 )' )
        else:
            if isRun2011Mu40eta2p1 or isRun2012:
                process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu40_eta2p1_v*",1,0 )' )
            elif isRun2011Mu40:
                process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu40_v*",1,0 )' )
    if runElectronAnalysis:
        switchOnTriggerMatching( process, triggerMatchers = [ 'electronTriggerMatchHLTElectrons' ], outputModule = '' )
        switchOnTriggerMatchEmbedding( process, triggerMatchers = [ 'electronTriggerMatchHLTElectrons' ], outputModule = '' )
        removeCleaningFromTriggerMatching( process, outputModule = '' )
        if isRun2012:
            process.electronTriggerMatchHLTElectrons.matchedCuts = cms.string( 'path( "HLT_*",1,0 )' )

## ================ ##
## Nominal Filter   ##
## ================ ##
process.hNu = cms.EDFilter(
    "HeavyNuFilter",
    muonTag      = cms.InputTag( 'selectedPatMuons' ),
    jetTag       = cms.InputTag( 'selectedPatJetsPFlow' ),
    electronTag  = cms.InputTag( 'selectedPatElectrons' ),
    minMu1pt     = cms.double(50.),
    minMu2pt     = cms.double(35.),
    minJetPt     = cms.double(25.),
    maxMuAbsEta  = cms.double(2.5),
    maxJetAbsEta = cms.double(2.6),
    minMuonJetdR = cms.double(0.3),

    muonTrackRelIsoLimit  = cms.double(0.1), # 10.0),
    maxVertexZsepCM       = cms.double(0.03),
    maxJetVZsepCM         = cms.double(0.1),

    isPFJets = cms.bool(True),

    heepVersion = cms.untracked.int32(heepVersion),

    mode = cms.int32(smode), #0 for Z+Jets, 1 for TTBar and other
    electronRho  = cms.InputTag( 'kt6PFJetsCentral','rho','RECO' )
    )

process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') )
process.out.outputCommands = cms.untracked.vstring("drop *", "keep *_*_*_RECO", "keep *_*_*_LHE", "keep *_*_*_SIM", "keep *_*_*_HLT", "keep *_*_*_PATSKIM",
                                                   "drop patTaus_*_*_PATSKIM", "drop recoPFCandidates_*_*_PATSKIM", "drop recoIsoDepositedmValueMap_*_*_PATSKIM")
process.p = cms.Path( process.AnalysisIntroSequence * process.hNu )
process.end = cms.EndPath(process.out)