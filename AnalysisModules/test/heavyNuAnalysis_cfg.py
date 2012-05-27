import FWCore.ParameterSet.Config as cms

from operator import isSequenceType
import os

#--- Analysis Mode ---#
#options are HNUMU, HNUE, TOP, QCD, CLO
#analysisMode = 'HNUE'

#--- Data/MC switch ---#
isMC=False
isData=not isMC

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
isRun2011Mu24       = False
isRun2011Mu40       = False
isRun2011Mu40eta2p1 = False
isRun2012           = True

#--- This flag must be set to ensure Monte Carlo gets the right corrections
isRun2011A          = False

#--- Flags for data taking era, which are set automatically ---#
#--- Possible options for dataEra: 20111 (2011A), 20112 (2011B), 20121 (2012) ---#
dataEra   = 20121 
pileupEra = 20121
if not isRun2012:
    if isRun2011A:
        dataEra   = 20111
        pileupEra = 20111
        if is44x: 
            pileupEra = 20113
    else:
        dataEra   = 20112
        pileupEra = 20112
        if is44x: 
            pileupEra = 20114


#--- Flags for nominal studies ---#
runMuonAnalysis     = True
runElectronAnalysis = True
systematics    = False
tagandprobe    = False
#doTriggerStudy = False
doTriggerStudy = True

#--- HEEP ID for electrons ---#
#--- Recognized values: 40 (2012), 31 or 32 (2011) ---#
heepVersion = 40

#--- Flags for Top studies ---#
topStudy      = True
studyTopScale = False

#--- Flags for QCD studies ---#
qcdStudy  = False
doDijet   = False
doQuadJet = False
doClosure = False

#--- Should always be True ---#
isPileupMC = True
isPFJets   = True

process = cms.Process("PAT");

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)

# source
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring('/store/mc/Summer12/TTJets_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S7_START52_V5-v1/0000/46565F4D-3981-E111-8C87-002618943945.root')
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
    
# Corrections to isolation for electrons
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation            = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

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
       process.hnuTotalKinematicsFilter * process.patDefaultSequence * process.patTrackSequence * process.myRefitMuonSequence * process.kt6PFJetsForIsolation
   )
else:
   process.AnalysisIntroSequence = cms.Sequence(
       process.patDefaultSequence * process.patTrackSequence * process.myRefitMuonSequence * process.kt6PFJetsForIsolation
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


#--- Electrons trigger analysis ---#
process.dumpEvContent = cms.EDAnalyzer("EventContentAnalyzer")
#
#  the analysis comes in from here
#
import HeavyNu.AnalysisModules.heavyNuEleTriggerEff_cfi

if runElectronAnalysis:
    if doTriggerStudy:
        process.hNuEtriggerEff                     = HeavyNu.AnalysisModules.heavyNuEleTriggerEff_cfi.HeavyNuEleTriggerEff.clone()
        process.hNuEtriggerEff.plotFolderName      = cms.string('inclusive')
        process.hNuEtriggerEff.numOfflJets         = cms.int32(  0 )

        process.hNuEtriggerEffOneJ                 = HeavyNu.AnalysisModules.heavyNuEleTriggerEff_cfi.HeavyNuEleTriggerEff.clone()
        process.hNuEtriggerEffOneJ.plotFolderName  = cms.string('OneJet')
        process.hNuEtriggerEffOneJ.numOfflJets     = cms.int32(  1 )

        process.hNuEtriggerEffTwoJ                 = HeavyNu.AnalysisModules.heavyNuEleTriggerEff_cfi.HeavyNuEleTriggerEff.clone()
        process.hNuEtriggerEffTwoJ.plotFolderName  = cms.string('TwoJet')
        process.hNuEtriggerEffTwoJ.numOfflJets     = cms.int32(  2 )

        process.TriggerStudyElectronSequence = cms.Sequence( process.kt6PFJetsForIsolation
                                                             #* process.dumpEvContent
                                                             * process.hNuEtriggerEff
                                                             * process.hNuEtriggerEffOneJ
                                                             * process.hNuEtriggerEffTwoJ
                                                             )


#--- Output histgram file ---#
process.TFileService = cms.Service("TFileService",
       fileName = cms.string("analysis.root"),
)

## ================ ##
## Nominal Analysis ##
## ================ ##
process.load("HeavyNu.AnalysisModules.heavynuanalysis_cfi")
if isMCsignal:
    process.hNu.isSignal = cms.bool(True)
    process.load("HeavyNu.AnalysisModules.heavyNuGenFilter_cfi")
    process.hNuGenFilter.keepIds = cms.vint32(2,)

process.hNu.minMu2pt         = cms.double(40.)
if not isRun2012:
    process.hNu.minMu2pt     = cms.double(30.)
process.hNu.ZmassWinMinGeV   = cms.double(60.)
process.hNu.ZmassWinMaxGeV   = cms.double(120.)
process.hNu.isPFJets         = cms.bool(isPFJets)
process.hNu.studyMuSelectEff = cms.bool(tagandprobe)
process.hNu.oneTPcand        = cms.bool(True)
process.hNu.studyScaleFactor = cms.bool(False)
process.hNu.studyRatePerRun  = cms.bool(isData)
process.hNu.electronRho      = cms.InputTag( 'kt6PFJetsForIsolation','rho' )
#--- Values below zero disable the vertex requirement ---#
process.hNu.maxVertexZsepCM = cms.double(-1)
process.hNu.maxJetVZsepCM   = cms.double(-1)
#--- Pileup corrections ---#
process.hNu.pileupEra = cms.int32(pileupEra)
#options are HNUMU, HNUE, TOP, QCD, CLO
process.hNu.analysisMode = cms.untracked.string('HNUMU')
process.hNu.heepVersion = cms.untracked.int32(heepVersion)

#--- Muon ID corrections are taken from Dec 11, 2011 studies---#
process.hNu.applyMuIDEffcorr = cms.bool(isMC)

if isData:
    if runMuonAnalysis:
        process.hNu.muonTag     = cms.InputTag( 'selectedPatMuonsTriggerMatch' )
    else:
        process.hNu.muonTag     = cms.InputTag( 'selectedPatMuons' )
    process.hNu.trigMatchPset.trigEventTag    = cms.InputTag("patTriggerEvent")
    process.hNu.trigMatchPset.electronFilters = cms.vstring('hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter')
    process.hNu.muIDPset.eraForId             = cms.int32( dataEra )
    if isRun2011Mu24:
        process.hNu.trigMatchPset.muonTriggers     = cms.vstring( 'HLT_Mu24_v1','HLT_Mu24_v2' )
        process.hNu.trigMatchPset.electronTriggers = cms.vstring( '' )
        process.hNu.trigMatchPset.triggerPt = cms.double( 24. )
    else:
        process.hNu.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_v1','HLT_Mu40_v2','HLT_Mu40_v3','HLT_Mu40_v5','HLT_Mu40_eta2p1_v1','HLT_Mu40_eta2p1_v2','HLT_Mu40_eta2p1_v3','HLT_Mu40_eta2p1_v4','HLT_Mu40_eta2p1_v5','HLT_Mu40_eta2p1_v6','HLT_Mu40_eta2p1_v7','HLT_Mu40_eta2p1_v8','HLT_Mu40_eta2p1_v9' ) 
        process.hNu.trigMatchPset.electronTriggers = cms.vstring( 'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3','HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v4','HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v5','HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v6' )
        process.hNu.trigMatchPset.triggerPt = cms.double( 40. )
else:
    # turn on MC trigger simulation
    process.hNu.trigMatchPset.randomSeed = cms.int32(os.getpid())
    process.hNu.trigMatchPset.trigEra    = cms.int32( dataEra )
    # Parameters for muon ID corrections to MC
    process.hNu.muIDPset.eraForId        = cms.int32( dataEra )

if isPFJets:
    process.hNu.jetTag = cms.InputTag( 'selectedPatJetsPFlow' )

if Training:
    process.hNu.trainingFileName=cms.untracked.string("changeme_nntraining.txt")

#--- Necessary to sort out trigger complications ---#
process.hNuMu24       = process.hNu.clone()
process.hNuMu40       = process.hNu.clone() 
process.hNuMu40eta2p1 = process.hNu.clone()

process.hNuE          = process.hNu.clone(analysisMode = cms.untracked.string('HNUE'))

process.hNuMu24.trigMatchPset.triggerPt    = cms.double( 24. )
process.hNuMu24.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu24_v1','HLT_Mu24_v2' )
process.hNuMu24.trigMatchPset.randomSeed   = cms.int32( os.getpid() )
process.hNuMu40.trigMatchPset.triggerPt    = cms.double( 40. )
process.hNuMu40.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_v1','HLT_Mu40_v2','HLT_Mu40_v3','HLT_Mu40_v5' ) 
process.hNuMu40.trigMatchPset.randomSeed   = cms.int32( os.getpid() )
process.hNuMu40eta2p1.trigMatchPset.triggerPt    = cms.double( 40. )
process.hNuMu40eta2p1.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_eta2p1_v1','HLT_Mu40_eta2p1_v2','HLT_Mu40_eta2p1_v3','HLT_Mu40_eta2p1_v4','HLT_Mu40_eta2p1_v5','HLT_Mu40_eta2p1_v6','HLT_Mu40_eta2p1_v7','HLT_Mu40_eta2p1_v8','HLT_Mu40_eta2p1_v9' ) 
process.hNuMu40eta2p1.trigMatchPset.randomSeed   = cms.int32( os.getpid() )

process.hNuMu24jesHi = process.hNuMu24.clone( applyJECUsign = cms.int32(1) )
process.hNuMu24jesLo = process.hNuMu24.clone( applyJECUsign = cms.int32(-1) )
process.hNuMu40jesHi = process.hNuMu40.clone( applyJECUsign = cms.int32(1) )
process.hNuMu40jesLo = process.hNuMu40.clone( applyJECUsign = cms.int32(-1) )

process.hNuEjesHi = process.hNuE.clone( analysisMode = cms.untracked.string('HNUE'), applyJECUsign = cms.int32(1) )
process.hNuEjesLo = process.hNuE.clone( analysisMode = cms.untracked.string('HNUE'), applyJECUsign = cms.int32(-1) )

process.hNuMu24jerHi = process.hNuMu24.clone( applyJERsign = cms.int32(1) )
process.hNuMu24jerLo = process.hNuMu24.clone( applyJERsign = cms.int32(-1) )
process.hNuMu40jerHi = process.hNuMu40.clone( applyJERsign = cms.int32(1) )
process.hNuMu40jerLo = process.hNuMu40.clone( applyJERsign = cms.int32(-1) )

process.hNuMu24mesHi = process.hNuMu24.clone( applyMESfactor = cms.double(1.01) )
process.hNuMu24mesLo = process.hNuMu24.clone( applyMESfactor = cms.double(0.99) )
process.hNuMu40mesHi = process.hNuMu40.clone( applyMESfactor = cms.double(1.01) )
process.hNuMu40mesLo = process.hNuMu40.clone( applyMESfactor = cms.double(0.99) )

process.hNuMu24mer = process.hNuMu24.clone( checkMERUnc = cms.bool(True) )
process.hNuMu40mer = process.hNuMu40.clone( checkMERUnc = cms.bool(True) )

process.hNuMu24midHi = process.hNuMu24.clone( applyMuIDEffsign = cms.int32(1) )
process.hNuMu24midLo = process.hNuMu24.clone( applyMuIDEffsign = cms.int32(-1) )
process.hNuMu40midHi = process.hNuMu40.clone( applyMuIDEffsign = cms.int32(1) )
process.hNuMu40midLo = process.hNuMu40.clone( applyMuIDEffsign = cms.int32(-1) )
process.hNuMu24midHi.applyMuIDEffcorr = cms.bool( isMC )
process.hNuMu24midLo.applyMuIDEffcorr = cms.bool( isMC )
process.hNuMu40midHi.applyMuIDEffcorr = cms.bool( isMC )
process.hNuMu40midLo.applyMuIDEffcorr = cms.bool( isMC )

process.hNuMu24trigHi = process.hNuMu24.clone( applyTrigEffsign  = cms.int32(1) )
process.hNuMu24trigLo = process.hNuMu24.clone( applyTrigEffsign  = cms.int32(-1) )
process.hNuMu40trigHi = process.hNuMu40.clone( applyTrigEffsign  = cms.int32(1) )
process.hNuMu40trigLo = process.hNuMu40.clone( applyTrigEffsign  = cms.int32(-1) )

# Pileup uncertainty: +/- 8% on number of interactions leads to 0.4 in 2011A, 0.7 in 2011B
process.hNuMu24puHi = process.hNuMu24.clone( systPileupShift = cms.double(0.4) )
process.hNuMu24puLo = process.hNuMu24.clone( systPileupShift = cms.double(-0.4) )
process.hNuMu40puHi = process.hNuMu40.clone( systPileupShift = cms.double(0.7) )
process.hNuMu40puLo = process.hNuMu40.clone( systPileupShift = cms.double(-0.7) )

process.hNuEpuHi = process.hNuE.clone( analysisMode = cms.untracked.string('HNUE'), systPileupShift = cms.double(0.7) )
process.hNuEpuLo = process.hNuE.clone( analysisMode = cms.untracked.string('HNUE'), systPileupShift = cms.double(-0.7) )

## ============ ##
## QCD Analysis ##
## ============ ##
process.load("HeavyNu.AnalysisModules.heavynuqcd_cfi")

process.hNuQCD.isPFJets = cms.bool(isPFJets)
process.hNuQCD.getSurvivalRate = cms.bool( doDijet )    
process.hNuQCD.doClosureTest   = cms.bool( doClosure )    
process.hNuQCD.doQuadJetTest   = cms.bool( doQuadJet )    
#--- Results from 2/fb 2011 data --> Run 2011A only ---#
process.hNuQCD.reweightPtLow  = cms.vdouble( 30,40,50,60,80,100,200 )
process.hNuQCD.reweightPtHigh = cms.vdouble( 40,50,60,80,100,200,1000 )
process.hNuQCD.reweight2011A  = cms.vdouble( 0.059932,0.0630225,0.0682759,0.0770292,0.0943841,0.113484,0.170732 )
process.hNuQCD.reweight2011B  = cms.vdouble( 0.0681898,0.0683649,0.0738487,0.0793934,0.0959666,0.11178,0.196970 )
#--- Values below zero disable the vertex requirement ---#
process.hNuQCD.dimuonMaxVertexZsepCM = cms.double( -1.0 )
process.hNuQCD.maxJetVZsepCM         = cms.double( -1.0 )
#--- Pileup corrections: needed even for data! ---#
process.hNuQCD.pileupEra = cms.int32(pileupEra)

if isData:
    process.hNuQCD.muonTag                     = cms.InputTag( 'selectedPatMuonsTriggerMatch' )
    process.hNuQCD.trigMatchPset.trigEventTag  = cms.InputTag( 'patTriggerEvent' )
    process.hNuQCD.trigMatchPset.electronFilters  = cms.vstring('') 
    if isRun2011Mu24:
        process.hNuQCD.trigMatchPset.muonTriggers     = cms.vstring( 'HLT_Mu24_v1','HLT_Mu24_v2' )
        process.hNuQCD.trigMatchPset.triggerPt        = cms.double( 24. )
        process.hNuQCD.trigMatchPset.electronTriggers = cms.vstring( '' )
    else:
        process.hNuQCD.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_v1','HLT_Mu40_v2','HLT_Mu40_v3','HLT_Mu40_v5','HLT_Mu40_eta2p1_v1','HLT_Mu40_eta2p1_v2','HLT_Mu40_eta2p1_v3','HLT_Mu40_eta2p1_v4','HLT_Mu40_eta2p1_v5','HLT_Mu40_eta2p1_v6','HLT_Mu40_eta2p1_v7','HLT_Mu40_eta2p1_v8','HLT_Mu40_eta2p1_v9' ) 
        process.hNuQCD.trigMatchPset.triggerPt    = cms.double( 40. )
        process.hNuQCD.trigMatchPset.electronTriggers = cms.vstring( '' )
else:
    # turn on MC trigger simulation
    process.hNuQCD.trigMatchPset.randomSeed = cms.int32( os.getpid() )
    process.hNuQCD.trigMatchPset.trigEra    = cms.int32( dataEra )
    # Parameters for muon ID corrections
    process.hNuQCD.muIDPset.eraForId        = cms.int32( dataEra )

if isPFJets:
    process.hNuQCD.jetTag  = cms.InputTag( 'selectedPatJetsPFlow' )

process.hNuQCDMu24 = process.hNuQCD.clone()
process.hNuQCDMu40 = process.hNuQCD.clone() 

process.hNuQCDMu24.trigMatchPset.triggerPt  = cms.double( 24. )
process.hNuQCDMu24.trigMatchPset.randomSeed = cms.int32( os.getpid() )
process.hNuQCDMu40.trigMatchPset.triggerPt  = cms.double( 40. )
process.hNuQCDMu40.trigMatchPset.randomSeed = cms.int32( os.getpid() )

## ============ ##
## Top Analysis ##
## ============ ##
process.load("HeavyNu.AnalysisModules.heavynutopanalysis_cfi")
# process.hNuTop.minLep2pt        = cms.double(30.)
process.hNuTop.studyScaleFactor = cms.bool(studyTopScale)
process.hNuTop.heepVersion      = cms.int32(heepVersion)
process.hNuTop.electronRho      = cms.InputTag( 'kt6PFJetsForIsolation','rho' )

#--- Some corrections re-enabled ---#
process.hNuTop.applyMuIDEffcorr      = cms.bool(isMC)
process.hNuTop.applyMuIDEffsign      = cms.int32(0)
process.hNuTop.applyEleEScale        = cms.bool(isMC) 
process.hNuTop.applyEleIDweight      = cms.bool(isMC) 
process.hNuTop.pileupEra             = cms.int32(pileupEra)
process.hNuTop.trigMatchPset.trigEra = cms.int32(dataEra)
process.hNuTop.muIDPset.eraForId     = cms.int32(dataEra)
process.hNuTop.EBidWgt = cms.double( 1.000 ) 
process.hNuTop.EEidWgt = cms.double( 1.000 ) 
#--- Values below zero disable the vertex requirement ---#
process.hNuTop.maxVertexZsepCM     = cms.double(-1)
process.hNuTop.maxVertexJetVZsepCM = cms.double(-1)

if isPFJets:
    process.hNuTop.jetTag = cms.InputTag( 'selectedPatJetsPFlow' )

if isData:
    process.hNuTop.muonTag                       = cms.InputTag( 'selectedPatMuonsTriggerMatch' )
    process.hNuTop.trigMatchPset.trigEventTag    = cms.InputTag("patTriggerEvent")
    process.hNuTop.trigMatchPset.electronFilters = cms.vstring('')
    if isRun2011Mu24:
        process.hNuTop.trigMatchPset.muonTriggers     = cms.vstring( 'HLT_Mu24_v1','HLT_Mu24_v2' )
        process.hNuTop.trigMatchPset.triggerPt        = cms.double( 24. )
        process.hNuTop.trigMatchPset.electronTriggers = cms.vstring( '' )
    else:
        process.hNuTop.trigMatchPset.muonTriggers     = cms.vstring( 'HLT_Mu40_v1','HLT_Mu40_v2','HLT_Mu40_v3','HLT_Mu40_v5','HLT_Mu40_eta2p1_v1','HLT_Mu40_eta2p1_v2','HLT_Mu40_eta2p1_v3','HLT_Mu40_eta2p1_v4','HLT_Mu40_eta2p1_v5','HLT_Mu40_eta2p1_v6','HLT_Mu40_eta2p1_v7','HLT_Mu40_eta2p1_v8','HLT_Mu40_eta2p1_v9' ) 
        process.hNuTop.trigMatchPset.triggerPt        = cms.double( 40. )
        process.hNuTop.trigMatchPset.electronTriggers = cms.vstring( '' )
else:
    process.hNuTop.trigMatchPset.randomSeed=cms.int32( os.getpid() )
    
process.hNuTopMu24 = process.hNuTop.clone()
process.hNuTopMu40 = process.hNuTop.clone() 

process.hNuTopMu24.trigMatchPset.triggerPt  = cms.double( 24. )
process.hNuTopMu24.trigMatchPset.randomSeed = cms.int32( os.getpid() )
process.hNuTopMu40.trigMatchPset.triggerPt  = cms.double( 40. )
process.hNuTopMu40.trigMatchPset.randomSeed = cms.int32( os.getpid() )

process.hNuTopMu24midHi = process.hNuTopMu24.clone( applyMuIDEffsign = cms.int32(1) )
process.hNuTopMu24midLo = process.hNuTopMu24.clone( applyMuIDEffsign = cms.int32(-1) )
process.hNuTopMu40midHi = process.hNuTopMu40.clone( applyMuIDEffsign = cms.int32(1) )
process.hNuTopMu40midLo = process.hNuTopMu40.clone( applyMuIDEffsign = cms.int32(-1) )

process.hNuTopMu24trigHi = process.hNuTopMu24.clone( applyTrigEffsign  = cms.int32(1) )
process.hNuTopMu24trigLo = process.hNuTopMu24.clone( applyTrigEffsign  = cms.int32(-1) )
process.hNuTopMu40trigHi = process.hNuTopMu40.clone( applyTrigEffsign  = cms.int32(1) )
process.hNuTopMu40trigLo = process.hNuTopMu40.clone( applyTrigEffsign  = cms.int32(-1) )

#-------------#
#--- Paths ---#
#-------------#
if runMuonAnalysis:
    if qcdStudy and topStudy:
        process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('pNominal','pQCD','pTop') )
    if qcdStudy and not topStudy:
        process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('pNominal','pQCD') )
    if topStudy and not qcdStudy:
        process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('pNominal','pTop') )
    if not qcdStudy and not topStudy:
        process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('pNominal') )
else: 
    if qcdStudy and topStudy:
        process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('pQCD','pTop') )
    if qcdStudy and not topStudy:
        process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('pQCD') )
    if topStudy and not qcdStudy:
        process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('pTop') )

if runMuonAnalysis:
    if isMC:
        if isMCsignal:
            process.AnalysisIntroSequence += process.hNuGenFilter

        if isRun2011A:
            process.p24 = cms.Path( process.AnalysisIntroSequence + process.hNuMu24 )
        process.p40 = cms.Path( process.AnalysisIntroSequence + process.hNuMu40 ) 

        if systematics:
            if isRun2011A:
                process.p24jesHi  = cms.Path( process.AnalysisIntroSequence + process.hNuMu24jesHi )
                process.p24jesLo  = cms.Path( process.AnalysisIntroSequence + process.hNuMu24jesLo )
                process.p24jerHi  = cms.Path( process.AnalysisIntroSequence + process.hNuMu24jerHi )
                process.p24jerLo  = cms.Path( process.AnalysisIntroSequence + process.hNuMu24jerLo )
                process.p24mesHi  = cms.Path( process.AnalysisIntroSequence + process.hNuMu24mesHi )
                process.p24mesLo  = cms.Path( process.AnalysisIntroSequence + process.hNuMu24mesLo )
                process.p24mer    = cms.Path( process.AnalysisIntroSequence + process.hNuMu24mer )
                process.p24midHi  = cms.Path( process.AnalysisIntroSequence + process.hNuMu24midHi )
                process.p24midLo  = cms.Path( process.AnalysisIntroSequence + process.hNuMu24midLo )
                process.p24trigHi = cms.Path( process.AnalysisIntroSequence + process.hNuMu24trigHi )
                process.p24trigLo = cms.Path( process.AnalysisIntroSequence + process.hNuMu24trigLo )
                process.p24puHi   = cms.Path( process.AnalysisIntroSequence + process.hNuMu24puHi )
                process.p24puLo   = cms.Path( process.AnalysisIntroSequence + process.hNuMu24puLo )

            process.p40jesHi  = cms.Path( process.AnalysisIntroSequence + process.hNuMu40jesHi )
            process.p40jesLo  = cms.Path( process.AnalysisIntroSequence + process.hNuMu40jesLo )
            process.p40jerHi  = cms.Path( process.AnalysisIntroSequence + process.hNuMu40jerHi )
            process.p40jerLo  = cms.Path( process.AnalysisIntroSequence + process.hNuMu40jerLo )
            process.p40mesHi  = cms.Path( process.AnalysisIntroSequence + process.hNuMu40mesHi )
            process.p40mesLo  = cms.Path( process.AnalysisIntroSequence + process.hNuMu40mesLo )
            process.p40mer    = cms.Path( process.AnalysisIntroSequence + process.hNuMu40mer )
            process.p40midHi  = cms.Path( process.AnalysisIntroSequence + process.hNuMu40midHi )
            process.p40midLo  = cms.Path( process.AnalysisIntroSequence + process.hNuMu40midLo )
            process.p40trigHi = cms.Path( process.AnalysisIntroSequence + process.hNuMu40trigHi )
            process.p40trigLo = cms.Path( process.AnalysisIntroSequence + process.hNuMu40trigLo )
            process.p40puHi   = cms.Path( process.AnalysisIntroSequence + process.hNuMu40puHi )
            process.p40puLo   = cms.Path( process.AnalysisIntroSequence + process.hNuMu40puLo )
    else:
        process.pNominal = cms.Path( process.AnalysisIntroSequence + process.hNu )
        if doTriggerStudy:
            if isRun2011Mu24:
                process.p24 = cms.Path( process.AnalysisIntroSequence + process.hNuMu24 )
            else:
                if isRun2011Mu40eta2p1 or isRun2012:
                    process.p40eta2p1 = cms.Path( process.AnalysisIntroSequence + process.hNuMu40eta2p1 )
                elif isRun2011Mu40:
                    process.p40       = cms.Path( process.AnalysisIntroSequence + process.hNuMu40 )

if runElectronAnalysis:
   if isMC:
      if isMCsignal:
         process.AnalysisIntroSequence += process.hNuGenFilter

   process.pE = cms.Path(process.AnalysisIntroSequence + process.hNuE)
   if systematics:
      process.pEjesHi = cms.Path(process.AnalysisIntroSequence + process.hNuEjesHi);
      process.pEjesLo = cms.Path(process.AnalysisIntroSequence + process.hNuEjesLo);
      process.pEpuHi = cms.Path(process.AnalysisIntroSequence + process.hNuEpuHi);
      process.pEpuLo = cms.Path(process.AnalysisIntroSequence + process.hNuEpuLo);
      
   if doTriggerStudy:
       process.pENominal = cms.Path(  process.TriggerStudyElectronSequence + process.AnalysisIntroSequence + process.hNuE  )
   else:
       process.pENominal = cms.Path( process.AnalysisIntroSequence + process.hNuE )
          
if qcdStudy:
    if isMC:
        if isRun2011A:
            process.pQCD24 = cms.Path( process.AnalysisIntroSequence + process.hNuQCDMu24 )
        process.pQCD40 = cms.Path( process.AnalysisIntroSequence + process.hNuQCDMu40 ) 
    else:
        process.pQCD = cms.Path( process.AnalysisIntroSequence + process.hNuQCD )
        
if topStudy:
    if isMC:
        if isRun2011A:
            process.pTop24 = cms.Path( process.AnalysisIntroSequence + process.hNuTopMu24 )
        process.pTop40 = cms.Path( process.AnalysisIntroSequence + process.hNuTopMu40 ) 
        if systematics:
            if isRun2011A:
                process.pTop24midHi = cms.Path( process.AnalysisIntroSequence + process.hNuTopMu24midHi )
                process.pTop24midLo = cms.Path( process.AnalysisIntroSequence + process.hNuTopMu24midLo )
                process.pTop24trigHi = cms.Path( process.AnalysisIntroSequence +  process.hNuTopMu24trigHi)
                process.pTop24trigLo = cms.Path( process.AnalysisIntroSequence +  process.hNuTopMu24trigLo)

            process.pTop40midHi = cms.Path( process.AnalysisIntroSequence + process.hNuTopMu40midHi )
            process.pTop40midLo = cms.Path( process.AnalysisIntroSequence + process.hNuTopMu40midLo )
            process.pTop40trigHi = cms.Path( process.AnalysisIntroSequence +  process.hNuTopMu40trigHi)
            process.pTop40trigLo = cms.Path( process.AnalysisIntroSequence +  process.hNuTopMu40trigLo)
    else:
        process.pTop = cms.Path( process.AnalysisIntroSequence + process.hNuTop )

