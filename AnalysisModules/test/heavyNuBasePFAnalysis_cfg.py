import FWCore.ParameterSet.Config as cms

import os

#--- Data/MC switch ---#
isMC=False
isData=not isMC

#--- Signal MC flags ---#
isMCsignal=False
Training=False

#--- Trigger-based luminosity flags ---#
#--- only one should be True        ---#
#--- LoLumi     --> HLT_Mu24        ---#
#--- HiLumi     --> HLT_Mu40        ---#
#--- VeryHiLumi --> HLT_Mu40_eta2p1 ---#
isRun2011LoLumi     = False
isRun2011HiLumi     = True
isRun2011VeryHiLumi = False

#--- Flags for data taking era ---#
isRun2011A = True
muIdYear   = 20110
pileupEra  = 20113
#--- Placeholder for 2011B variables
if not isRun2011A:
    muIdYear  = 20111
    pileupEra = 20113


runAnalysis = True

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
    fileNames=cms.untracked.vstring('file:/hdfs/cms/user/dahmes/storage/wr2011/MuSkim/ReReco/may10/may10rerecoMuSkim_022.root')
)

if isData:
    if isRun2011LoLumi:
        print "===========> Flag is SET for 2011 LOW luminosity data <============"
        from HeavyNu.AnalysisModules.goodLumiList_160404_163869_may10rereco_Mu24_cfi import lumisToProcess
    else:
        if isRun2011VeryHiLumi:
            print "===========> Flag is SET for 2011 HIGH luminosity data <============"
            from HeavyNu.AnalysisModules.goodLumiList_173236_180252_Mu40eta2p1_cfi import lumisToProcess
        else:
            print "===========> Flag is SET for 2011 MEDIUM luminosity data <============"
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
        process.GlobalTag.globaltag=cms.string('START42_V13::All')
        print "=============> isPileupMC flag is SET <================"
    # else:
    #     print "========> Fall10 MC with Spring10 JEC applied <========"
    #     process.GlobalTag.globaltag = cms.string('START38_V14::All')
else:
    print "===============> Running on DATA <===================="
    process.GlobalTag.globaltag = cms.string('GR_R_42_V20::All')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')

################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

## pat sequences to be loaded:
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("RecoMuon.MuonIdentification.refitMuons_cfi")
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
    # Needed for AOD PF jets
    ##-------------------- Import the JEC services -----------------------
    process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
    process.kt6PFJets = kt4PFJets.clone(
        rParam = cms.double(0.6),
        src = cms.InputTag('pfNoElectron'),
        doAreaFastjet = cms.bool(True),
        doRhoFastjet = cms.bool(True)
    )

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
    # process.outpath  = cms.EndPath(process.out)
    from PhysicsTools.PatAlgos.tools.coreTools import *
    if isPFJets:
        removeMCMatchingPF2PAT( process, '' )
    removeMCMatching(process, ['All'], outputInProcess = False)
        

#--- Calo Jet Energy Corrections: No longer used ---#
process.patJetCorrFactors.useRho = cms.bool(True)
#Corrections for Calo Jets: no longer used
if isMC:
    switchJetCollection( process,
                         jetCollection=cms.InputTag('ak5PFJets',"","RECO"),
                         jetCorrLabel=('AK5PF', ['L1FastJet','L2Relative','L3Absolute']))
else:
    switchJetCollection( process,
                         jetCollection=cms.InputTag('ak5PFJets', "", "RECO"),
                         jetCorrLabel=('AK5PF', ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']))

# Compute the mean pt per unit area (rho) from the PFchs inputs
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsWithRho = kt4PFJets.clone(
        rParam = cms.double(0.6),
        src = cms.InputTag('pfNoElectron'),
        doAreaFastjet = cms.bool(True),
        doRhoFastjet = cms.bool(True)
    )
process.patJetCorrFactors.rho = cms.InputTag("kt6PFJetsWithRho", "rho")

#now apparently we reconstruct PF2PAT just for the JEC
process.load("CommonTools/ParticleFlow/pfElectrons_cff")
process.load("CommonTools/ParticleFlow/pfMuons_cff")
process.load("CommonTools/ParticleFlow/pfNoPileUp_cff")
process.load("CommonTools/ParticleFlow/TopProjectors/pfNoMuon_cfi")
process.load("CommonTools/ParticleFlow/TopProjectors/pfNoElectron_cfi")

#Add new kt4PFJets collection with fastjet area to path
#process.patJetMETCorrections = cms.Sequence(cms.Sequence(process.kt6PFJets)*process.patJetCorrFactors)
getattr(process,"patDefaultSequence").replace(
        getattr(process,"patJetCorrections"),
        process.pfNoPileUpSequence*process.pfMuonSequence*process.pfNoMuon*process.pfElectronSequence*
         process.pfNoElectron*process.kt6PFJetsWithRho*getattr(process,"patJetCorrections")
    )

#process.patDefaultSequenceModified = cms.Sequence(
#        process.goodOfflinePrimaryVertices*
#        getattr(process,"patPF2PATSequence"+postfix)
#    )

## --------------------- ##
## Define the basic path ##
## --------------------- ##
process.AnalysisIntroSequence = cms.Sequence(
    process.patDefaultSequence * process.refitMuons
)
if isPFJets:
    process.AnalysisIntroSequence += process.modifiedPF2PATSequence
    # Needed for jet corrections with AOD PF jets 
    process.AnalysisIntroSequence += process.kt6PFJets

process.p = cms.Path(
    process.AnalysisIntroSequence
)


## ============================== ##
## Python tools --> Order matters ##
## ============================== ##

from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
removeCleaning( process, False )

# Special change for saving good products in data
if isData:
    process.out.outputCommands = cms.untracked.vstring("keep *")

#--- Trigger matching ---#
process.load("HeavyNu.AnalysisModules.hnutrigmatch_cfi")
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

# process.hNu.minMu2pt         = cms.double(30.)
process.hNu.isPFJets         = cms.bool(isPFJets)
process.hNu.studyMuSelectEff = cms.bool(False)
process.hNu.studyScaleFactor = cms.bool(False)
process.hNu.studyRatePerRun  = cms.bool(isData)
#--- Values below zero disable the vertex requirement ---#
process.hNu.maxVertexZsepCM = cms.double(-1)
process.hNu.maxJetVZsepCM   = cms.double(-1)
#--- Pileup corrections ---#
process.hNu.pileupEra = cms.int32(pileupEra)

#--- Muon ID corrections are available as of 27 Oct, 2011 ---#
process.hNu.applyMuIDEffcorr = cms.bool(isMC)

if isData:
    process.hNu.muonTag                    = cms.InputTag( 'selectedPatMuonsTriggerMatch' )
    process.hNu.trigMatchPset.trigEventTag = cms.InputTag("patTriggerEvent")
    process.hNu.trigMatchPset.muonMatch    = cms.string('muonTriggerMatchHLTMuons')
    if isRun2011LoLumi:
        process.hNu.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu24_v1','HLT_Mu24_v2' )
        process.hNu.trigMatchPset.triggerPt = cms.double( 24. )
    else:
        process.hNu.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_v1','HLT_Mu40_v2','HLT_Mu40_v3','HLT_Mu40_v5','HLT_Mu40_eta2p1_v1' ) 
        process.hNu.trigMatchPset.triggerPt = cms.double( 40. )
else:
    # turn on MC trigger simulation
    process.hNu.trigMatchPset.randomSeed = cms.int32(os.getpid())
    # Parameters for muon ID corrections to MC
    process.hNu.muIDPset.eraForId        = cms.int32( muIdYear )

if isPFJets:
    process.hNu.jetTag = cms.InputTag( 'selectedPatJetsPFlow' )

if Training:
    process.hNu.trainingFileName=cms.untracked.string("changeme_nntraining.txt")

if runAnalysis:
    if isMC:
        if isMCsignal:
            process.AnalysisIntroSequence += process.hNuGenFilter
    else:
        process.pNominal = cms.Path( process.AnalysisIntroSequence + process.hNu )

#terrible analysis with base jet collection
process.hNuBasePFJets = process.hNu.clone()
process.hNuBasePFJets.jetTag = cms.InputTag( 'selectedPatJets' )

if not isMC:
      process.pBasePFJets = cms.Path( process.AnalysisIntroSequence + process.hNuBasePFJets )

process.hNuFNAL = cms.EDFilter( "HeavyNuFNAL",
    trigMatchPset = cms.PSet(
        trigEventTag = cms.InputTag( "" ),
        muonTriggers = cms.vstring( '' ),
        firstRun     = cms.vint32( 0 ),
        lastRun      = cms.vint32( 999999 ),
        muonMatch    = cms.string( '' ),
        triggerPt    = cms.double( 40. ),
        randomSeed   = cms.int32( 0 ),  # for MC
        year         = cms.int32( 2011 ) # for MC
    ),
    DoLog        = cms.bool( False ),
    muonTag      = cms.InputTag( 'selectedPatMuons' ),
    jetTag       = cms.InputTag( 'selectedPatJets' ),
    metTag       = cms.InputTag( 'patMETs' ),
    BtagName     = cms.string( 'trackCountingHighEffBJetTags' ),
    minBtagDiscr = cms.double(2.0), 
    minMu1pt     = cms.double(25.),
    minMu2pt     = cms.double(20.),
    minJet1Pt    = cms.double(120),
    minJet2Pt    = cms.double(25),
    maxMuAbsEta  = cms.double(2.4),
    maxJetAbsEta = cms.double(3.0),
    minMuonJetdR = cms.double(0.5),
    muonTrackRelIsoLimit  = cms.double(999999.0), 
    maxVertexZsepCM       = cms.double(-1.),
    maxJetVZsepCM         = cms.double(-1.),

    jetsFromReco  = cms.bool( True ),
    jetCorrString = cms.string( "ak5PFL1FastL2L3Residual" ),

    useTrackerPt = cms.bool(True),
    isPFJets = cms.bool(False)
    )

#--- Runs baseline Fermilab configuration, including AOD jets ---#
process.pFNAL = cms.Path( process.AnalysisIntroSequence + process.hNuFNAL )

#--- Each subsequent selection adds cuts to BASE case ---#

#--- Switch to using cocktail muon pT, PAT jets from AOD input ---#
process.hNuFNALCocktail = process.hNuFNAL.clone( useTrackerPt = cms.bool(False), jetsFromReco = cms.bool(False) )
process.pFNALCocktail   = cms.Path( process.AnalysisIntroSequence + process.hNuFNALCocktail )

#--- Add jets selected using PF2PAT ---#
process.hNuFNALPF2PATCocktail = process.hNuFNALCocktail.clone( jetTag = cms.InputTag('selectedPatJetsPFlow') )
process.pFNALPF2PATCocktail   = cms.Path( process.AnalysisIntroSequence + process.hNuFNALPF2PATCocktail )

#--- Now try with jets within tracker ---#
process.hNuFNALeta2p5 = process.hNuFNALPF2PATCocktail.clone( maxJetAbsEta = cms.double(2.5) )
process.pFNALeta2p5   = cms.Path( process.AnalysisIntroSequence + process.hNuFNALeta2p5 )

#--- And finally with muon isolation ---#
process.hNuFNALmuIso = process.hNuFNALeta2p5.clone( muonTrackRelIsoLimit = cms.double(0.10) )
process.pFNALmuIso   = cms.Path( process.AnalysisIntroSequence + process.hNuFNALmuIso )

