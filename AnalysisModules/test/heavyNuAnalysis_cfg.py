import FWCore.ParameterSet.Config as cms

from operator import isSequenceType
import os

#--- Data/MC switch ---#
isMC=False
isData=not isMC

#--- Specify CMSSW release (53, 52, 44, 42, 37) ---#
cmsswRelease = 53
isRereco = True

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
#--- Special flag to use the *full* 5/fb for ICHEP ---#
#--- Default (i.e. False) is to use the 3.6/fb specified by EXO HN as acceptable ---#
isRun2012topoff     = False 

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
        if cmsswRelease == 44: 
            pileupEra = 20113
    else:
        dataEra   = 20112
        pileupEra = 20112
        if cmsswRelease == 44: 
            pileupEra = 20114

#--- Flags for nominal studies ---#
runMuonAnalysis     = True
runElectronAnalysis = True
systematics    = False
tagandprobe    = False
doTriggerStudy = False
addSlopeTrees  = True

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

#--- Tau Flag ---#
ishpsPFTau = False

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
                            fileNames=cms.untracked.vstring('file:/hdfs/cms/skim/mu/hNu_2012/jul13_2012A_mu35e35/jul13_2012A_mu35e35_015.root')
                            #/store/mc/Summer12_DR53X/WRToNuLeptonToLLJJ_MW-2900_MNu-1450_TuneZ2star_8TeV-pythia6-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/4005DF3C-ABEC-E111-BC0E-00215E21D570.root')
                            #file:/local/cms/user/pastika/heavyNuAnalysis_2012/2C1FBAB2-C1D4-E111-A89A-001E6739815B.root
                            #file:/home/ugrad/pastika/cms/HeavyNu/CMSSW_5_3_3_patch1/src/HeavyNu/AnalysisModules/heavynu_candevents.root')
                            #file:/hdfs/cms/skim/mu/hNu_2012/jul13_2012A_mu35e35/jul13_2012A_mu35e35_015.root
)

if isData:
    if isRun2012:
        if isRereco:
            print "===========> Flag is SET for Reprocessing 5/fb 2012 data <============"
            from HeavyNu.AnalysisModules.goodLumiList_2012_reReco_July13_cfi import lumisToProcess
        elif isRun2012topoff:
            print "===========> Flag is SET for 5/fb 2012 data <============"
            from HeavyNu.AnalysisModules.goodLumiList_2012_dynamic_topoff_ichep_cfi import lumisToProcess
        else:
            print "===========> Flag is SET for 3.6/fb 2012 data <============"
            from HeavyNu.AnalysisModules.goodLumiList_2012_dynamic_cfi import lumisToProcess
    else:
        if isRun2011Mu24:
            print "===========> Flag is SET for 2011 LOW luminosity data <============"
            if cmsswRelease == 44:
                from HeavyNu.AnalysisModules.goodLumiList_novRereco_160404_163869_mu24_cfi import lumisToProcess
            else:
                from HeavyNu.AnalysisModules.goodLumiList_160404_163869_may10rereco_Mu24_cfi import lumisToProcess
        else:
            if isRun2011Mu40eta2p1:
                print "===========> Flag is SET for 2011 HIGH luminosity data <============"
                if cmsswRelease == 44:
                    from HeavyNu.AnalysisModules.goodLumiList_novRereco_173236_180252_mu40_eta2p1_cfi import lumisToProcess
                else:
                    from HeavyNu.AnalysisModules.goodLumiList_173236_180252_Mu40eta2p1_cfi import lumisToProcess
            elif isRun2011Mu40:
                print "===========> Flag is SET for 2011 MEDIUM luminosity data <============"
                if cmsswRelease == 44:
                    from HeavyNu.AnalysisModules.goodLumiList_novRereco_165088_173198_mu40_cfi import lumisToProcess
                else:
                    from HeavyNu.AnalysisModules.goodLumiList_165088_173198_Mu40_cfi import lumisToProcess
                            
    process.source.lumisToProcess = lumisToProcess

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


## Load additional processes
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

## Global Tags:
if (isMC):
    print "=================> MC flag is SET <===================="
    if cmsswRelease == 42:
        process.GlobalTag.globaltag=cms.string('START44_V12::All')
    elif cmsswRelease == 44:
        process.GlobalTag.globaltag=cms.string('START44_V12::All')
    elif cmsswRelease == 52:
        process.GlobalTag.globaltag=cms.string('START52_V11::All')
    elif cmsswRelease == 53:
        process.GlobalTag.globaltag=cms.string('START53_V7A::All')
    else:
        print "INVALID CMSSW release id %(rid)i"%{"rid":cmsswRelease}
else:
    print "===============> Running on DATA <===================="
    if cmsswRelease == 42:
        process.GlobalTag.globaltag = cms.string('GR_R_42_V20::All')
    elif cmsswRelease == 44:
        process.GlobalTag.globaltag = cms.string('GR_R_44_V13::All')
    elif cmsswRelease == 52:
        process.GlobalTag.globaltag = cms.string('GR_R_52_V9::All')
    elif cmsswRelease == 53:
        if isRereco:
            process.GlobalTag.globaltag = cms.string('FT_53_V6_AN2::All')
        else:
            process.GlobalTag.globaltag = cms.string('GR_P_V40::All')
    else:
        print "INVALID CMSSW release id %(rid)i"%{"rid":cmsswRelease}

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

def usePF2PAT_WREdition(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix="", jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute']), pvCollection=cms.InputTag('offlinePrimaryVertices'), typeIMetCorrections=False, outputModules=['out']):

    # -------- CORE ---------------
    if runPF2PAT:
        process.load("CommonTools.ParticleFlow.PF2PAT_cff")
        #add Pf2PAT *before* cloning so that overlapping modules are cloned too
        #process.patDefaultSequence.replace( process.patCandidates, process.PF2PAT+process.patCandidates)
        process.patPF2PATSequence = cms.Sequence( process.PF2PAT + process.patDefaultSequence)
    else:
        process.patPF2PATSequence = cms.Sequence( process.patDefaultSequence )

    if not postfix == "":
        from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet
        cloneProcessingSnippet(process, process.patPF2PATSequence, postfix)
        #delete everything pat PF2PAT modules! if you want to test the postfixing for completeness
        #from PhysicsTools.PatAlgos.tools.helpers import listModules,listSequences
        #for module in listModules(process.patDefaultSequence):
        #    if not module.label() is None: process.__delattr__(module.label())
        #for sequence in listSequences(process.patDefaultSequence):
        #    if not sequence.label() is None: process.__delattr__(sequence.label())
        #del process.patDefaultSequence

    removeCleaning(process, postfix=postfix, outputModules=outputModules)

    # -------- OBJECTS ------------
    # Muons
    #adaptPFMuons(process,
    #             applyPostfix(process,"patMuons",postfix),
    #             postfix)

    # Electrons
    #adaptPFElectrons(process,
    #                 applyPostfix(process,"patElectrons",postfix),
    #                 postfix)

    # Photons
    print "Temporarily switching off photons completely"

    removeSpecificPATObjects(process,names=['Photons', 'Taus', "Muons", "Electrons"],outputModules=outputModules,postfix=postfix)
    removeIfInSequence(process,"patPhotonIsolation","patDefaultSequence",postfix)

    # Jets
    if runOnMC :
       if ishpsPFTau:
         switchToPFJets( process, cms.InputTag('pfNoTau'+postfix), jetAlgo, postfix=postfix,
                        jetCorrections=jetCorrections, type1=typeIMetCorrections, outputModules=outputModules )
       else:
         switchToPFJets( process, cms.InputTag('pfJets'+postfix), jetAlgo, postfix=postfix,
                           jetCorrections=jetCorrections, type1=typeIMetCorrections, outputModules=outputModules )
         applyPostfix(process,"patDefaultSequence",postfix).replace(
            applyPostfix(process,"patJetGenJetMatch",postfix),
            getattr(process,"genForPF2PATSequence") *
            applyPostfix(process,"patJetGenJetMatch",postfix)
            )
    else :
        if not 'L2L3Residual' in jetCorrections[1]:
            print '#################################################'
            print 'WARNING! Not using L2L3Residual but this is data.'
            print 'If this is okay with you, disregard this message.'
            print '#################################################'
        if ishpsPFTau:
            switchToPFJets( process, cms.InputTag('pfNoTau'+postfix), jetAlgo, postfix=postfix,
                            jetCorrections=jetCorrections, type1=typeIMetCorrections, outputModules=outputModules )
        else:
            switchToPFJets( process, cms.InputTag('pfJets'+postfix), jetAlgo, postfix=postfix,
                            jetCorrections=jetCorrections, type1=typeIMetCorrections, outputModules=outputModules )

    # Taus
    #adaptPFTaus( process, tauType='shrinkingConePFTau', postfix=postfix )
    #adaptPFTaus( process, tauType='fixedConePFTau', postfix=postfix )
    if ishpsPFTau:
        adaptPFTaus( process, tauType='hpsPFTau', postfix=postfix )

    # MET
    switchToPFMET(process, cms.InputTag('pfMET'+postfix), type1=typeIMetCorrections, postfix=postfix)

    # Unmasked PFCandidates
    addPFCandidates(process,cms.InputTag('pfNoJet'+postfix),patLabel='PFParticles'+postfix,cut="",postfix=postfix)

    # adapt primary vertex collection
    adaptPVs(process, pvCollection=pvCollection, postfix=postfix)

    if runOnMC:
        process.load("CommonTools.ParticleFlow.genForPF2PAT_cff")
        getattr(process, "patDefaultSequence"+postfix).replace(
            applyPostfix(process,"patCandidates",postfix),
            process.genForPF2PATSequence+applyPostfix(process,"patCandidates",postfix)
            )
    else:
        removeMCMatchingPF2PAT(process,postfix=postfix,outputModules=outputModules)

    print "Done: PF2PAT interfaced to PAT, postfix=", postfix


#--- Output module: 
#--- Must be defined before PAT python tools will work
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('heavynu_candevents.root'),
                               # save only events passing the full path
                               SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('p')),
                               outputCommands = cms.untracked.vstring("keep *")
                               )
if isPFJets:
    from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
    process.goodOfflinePrimaryVertices = cms.EDFilter(
        "PrimaryVertexObjectFilter",
        filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
        src=cms.InputTag('offlinePrimaryVertices')
    )

    postfix = "PFlow"
    if isMC:
        usePF2PAT_WREdition(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=isMC, postfix=postfix,
                  jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute']),
                  pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))
    else:
        if isRun2012:
            usePF2PAT_WREdition(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=isMC, postfix=postfix,
                      jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']),
                  pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))
        else:
            usePF2PAT_WREdition(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=isMC, postfix=postfix,
                      jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']),
                  pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))
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

    process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
    process.load("RecoBTag.Configuration.RecoBTag_cff")
    process.modifiedPF2PATSequence = cms.Sequence(    
        process.goodOfflinePrimaryVertices*process.inclusiveVertexing * process.btagging *
        getattr(process,"patPF2PATSequence"+postfix)
    )

    # Remove unneeded tau sequences, some of which are very time consuming
    if not ishpsPFTau:
       process.modifiedPF2PATSequence.remove(getattr(process,"tauIsoDepositPFCandidates"+postfix) )
       process.modifiedPF2PATSequence.remove(getattr(process,"tauIsoDepositPFChargedHadrons"+postfix) )
       process.modifiedPF2PATSequence.remove(getattr(process,"tauIsoDepositPFNeutralHadrons"+postfix) )
       process.modifiedPF2PATSequence.remove(getattr(process,"tauIsoDepositPFGammas"+postfix) )
       process.modifiedPF2PATSequence.remove(getattr(process,"tauMatch"+postfix) )
       process.modifiedPF2PATSequence.remove(getattr(process,"tauGenJets"+postfix) )
       process.modifiedPF2PATSequence.remove(getattr(process,"tauGenJetsSelectorAllHadrons"+postfix) )
       process.modifiedPF2PATSequence.remove(getattr(process,"tauGenJetMatch"+postfix) )

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

## ---- ##
## Btag ##
## ---- ##

## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagJetProbabilityCalibration?redirectedfrom=CMS.SWGuideBTagJetProbabilityCalibration#Calibration_in_52x_and_53x_Data

if cmsswRelease == 52:
    process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
          tag = cms.string("TrackProbabilityCalibration_2D_2012DataTOT_v1_offline"),
          connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
          tag = cms.string("TrackProbabilityCalibration_3D_2012DataTOT_v1_offline"),
          connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
    )


## --------------------- ##
## Define the basic path ##
## --------------------- ##

#--- Beam background removal ---#
process.scrapingFilter      = cms.EDFilter("FilterOutScraping",
                                           applyfilter = cms.untracked.bool(True),
                                           debugOn = cms.untracked.bool(False),
                                           numtrack = cms.untracked.uint32(10),
                                           thresh = cms.untracked.double(0.25)
                                           )
#--- Primary vertex requirement ---#
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24), 
                                           maxd0 = cms.double(2) 
                                           )
#--- HB/HE event-level noise filter ---#
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi') 

process.eventFilters = cms.Sequence( process.scrapingFilter + process.primaryVertexFilter + process.HBHENoiseFilter ) 

if isMC:
   # Gen Level Energy balance filter to fix Pythia6 lhe interface bug
   process.load("HeavyNu.AnalysisModules.hnuTotalKinematicsFilter_cfi")
   process.AnalysisIntroSequence = cms.Sequence(
       process.hnuTotalKinematicsFilter * process.eventFilters * process.patDefaultSequence * process.patTrackSequence * process.myRefitMuonSequence * process.kt6PFJetsForIsolation
   )
else:
   process.AnalysisIntroSequence = cms.Sequence(
       process.eventFilters * process.patDefaultSequence * process.patTrackSequence * process.myRefitMuonSequence * process.kt6PFJetsForIsolation
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
    if cmsswRelease == 42:
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
from PhysicsTools.PatAlgos.tools.coreTools import removeSpecificPATObjects
if (cmsswRelease == 42):
    removeCleaning( process, False )
elif (cmsswRelease == 44):
    removeCleaning( process, outputModules = [] )
else:
    removeCleaning( process, outputModules = [] )
    removeSpecificPATObjects(process, names = ['Jets','Taus','METs'], outputModules = [])

# Special change for saving good products in data
if isData:
    process.out.outputCommands = cms.untracked.vstring("keep *")

#--- Trigger matching ---#
process.load("HeavyNu.AnalysisModules.hnutrigmatch_cfi")
from PhysicsTools.PatAlgos.tools.trigTools import *
if isData:
    if runMuonAnalysis or topStudy:
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

if runElectronAnalysis:
    if doTriggerStudy and isData:
        process.load("HeavyNu.AnalysisModules.heavyNuEleTriggerEff_cff")

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
process.hNu.oneTPcand        = cms.bool(False)
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
process.hNu.heepVersion  = cms.untracked.int32(heepVersion)
process.hNu.addSlopeTree = cms.untracked.bool(addSlopeTrees)
process.hNu.nFakeLeptons = cms.untracked.int32(0)

#--- Muon ID corrections are taken from June 22, 2012 studies---#
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

process.hNuE               = process.hNu.clone(analysisMode = cms.untracked.string('HNUE'))
process.hNuE.correctEscale = cms.bool(isMC)

#--- Electron ID corrections are taken from June 22, 2012 studies---#
process.hNu.applyMuIDEffcorr = cms.bool(isMC)

process.hNuMu24.trigMatchPset.triggerPt    = cms.double( 24. )
process.hNuMu24.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu24_v1','HLT_Mu24_v2' )
process.hNuMu24.trigMatchPset.randomSeed   = cms.int32( os.getpid() )
process.hNuMu40.trigMatchPset.triggerPt    = cms.double( 40. )
process.hNuMu40.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_v1','HLT_Mu40_v2','HLT_Mu40_v3','HLT_Mu40_v5' ) 
process.hNuMu40.trigMatchPset.randomSeed   = cms.int32( os.getpid() )
process.hNuMu40eta2p1.trigMatchPset.triggerPt    = cms.double( 40. )
process.hNuMu40eta2p1.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_eta2p1_v1','HLT_Mu40_eta2p1_v2','HLT_Mu40_eta2p1_v3','HLT_Mu40_eta2p1_v4','HLT_Mu40_eta2p1_v5','HLT_Mu40_eta2p1_v6','HLT_Mu40_eta2p1_v7','HLT_Mu40_eta2p1_v8','HLT_Mu40_eta2p1_v9' ) 
process.hNuMu40eta2p1.trigMatchPset.randomSeed   = cms.int32( os.getpid() )

process.hNuMu24jesHi = process.hNuMu24.clone( studyMuSelectEff = cms.bool(False), applyJECUsign = cms.int32(1) )
process.hNuMu24jesLo = process.hNuMu24.clone( studyMuSelectEff = cms.bool(False), applyJECUsign = cms.int32(-1) )
process.hNuMu40jesHi = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyJECUsign = cms.int32(1) )
process.hNuMu40jesLo = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyJECUsign = cms.int32(-1) )

process.hNuMu24jerHi = process.hNuMu24.clone( studyMuSelectEff = cms.bool(False), applyJERsign = cms.int32(1) )
process.hNuMu24jerLo = process.hNuMu24.clone( studyMuSelectEff = cms.bool(False), applyJERsign = cms.int32(-1) )
process.hNuMu40jerHi = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyJERsign = cms.int32(1) )
process.hNuMu40jerLo = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyJERsign = cms.int32(-1) )

process.hNuEjesHi  = process.hNuE.clone( studyMuSelectEff = cms.bool(False), applyJECUsign = cms.int32(1) )
process.hNuEjesLo  = process.hNuE.clone( studyMuSelectEff = cms.bool(False), applyJECUsign = cms.int32(-1) )
process.hNuEjerHi  = process.hNuE.clone( studyMuSelectEff = cms.bool(False), applyJERsign = cms.int32(1) )
process.hNuEjerLo  = process.hNuE.clone( studyMuSelectEff = cms.bool(False), applyJERsign = cms.int32(-1) )
process.hNuEescale = process.hNuE.clone( studyMuSelectEff = cms.bool(False), correctEscale = cms.bool(False) )
process.hNuEidHi   = process.hNuE.clone( studyMuSelectEff = cms.bool(False), applyMuIDEffsign = cms.int32(1) )
process.hNuEidLo   = process.hNuE.clone( studyMuSelectEff = cms.bool(False), applyMuIDEffsign = cms.int32(-1) )
process.hNuEidHi.applyMuIDEffcorr = cms.bool( isMC )
process.hNuEidLo.applyMuIDEffcorr = cms.bool( isMC )
process.hNuEtrigHi = process.hNuE.clone( studyMuSelectEff = cms.bool(False), applyTrigEffsign  = cms.int32(1) )
process.hNuEtrigLo = process.hNuE.clone( studyMuSelectEff = cms.bool(False), applyTrigEffsign  = cms.int32(-1) )

# process.hNuMu24mesHi = process.hNuMu24.clone( applyMESfactor = cms.double(1.01) )
# process.hNuMu24mesLo = process.hNuMu24.clone( applyMESfactor = cms.double(0.99) )
# process.hNuMu40mesHi = process.hNuMu40.clone( applyMESfactor = cms.double(1.01) )
# process.hNuMu40mesLo = process.hNuMu40.clone( applyMESfactor = cms.double(0.99) )

process.hNuMu24mer = process.hNuMu24.clone( studyMuSelectEff = cms.bool(False), checkMERUnc = cms.bool(True) )
process.hNuMu40mer = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), checkMERUnc = cms.bool(True) )

process.hNuMu24midHi = process.hNuMu24.clone( studyMuSelectEff = cms.bool(False), applyMuIDEffsign = cms.int32(1) )
process.hNuMu24midLo = process.hNuMu24.clone( studyMuSelectEff = cms.bool(False), applyMuIDEffsign = cms.int32(-1) )
process.hNuMu40midHi = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyMuIDEffsign = cms.int32(1) )
process.hNuMu40midLo = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyMuIDEffsign = cms.int32(-1) )
process.hNuMu24midHi.applyMuIDEffcorr = cms.bool( isMC )
process.hNuMu24midLo.applyMuIDEffcorr = cms.bool( isMC )
process.hNuMu40midHi.applyMuIDEffcorr = cms.bool( isMC )
process.hNuMu40midLo.applyMuIDEffcorr = cms.bool( isMC )

process.hNuMu24trigHi = process.hNuMu24.clone( studyMuSelectEff = cms.bool(False), applyTrigEffsign  = cms.int32(1) )
process.hNuMu24trigLo = process.hNuMu24.clone( studyMuSelectEff = cms.bool(False), applyTrigEffsign  = cms.int32(-1) )
process.hNuMu40trigHi = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyTrigEffsign  = cms.int32(1) )
process.hNuMu40trigLo = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyTrigEffsign  = cms.int32(-1) )

# Pileup uncertainty: +/- 8% on number of interactions leads to 0.4 in 2011A, 0.7 in 2011B
#                     +/- 5% on number of interactions (16.85, 12/06/09) leads to 0.84 in 2012AB
if isRun2012:
    process.hNuMu40puHi = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), systPileupShift = cms.double(0.84) )
    process.hNuMu40puLo = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), systPileupShift = cms.double(-0.84) )
else:
    process.hNuMu24puHi = process.hNuMu24.clone( studyMuSelectEff = cms.bool(False), systPileupShift = cms.double(0.4) )
    process.hNuMu24puLo = process.hNuMu24.clone( studyMuSelectEff = cms.bool(False), systPileupShift = cms.double(-0.4) )
    process.hNuMu40puHi = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), systPileupShift = cms.double(0.7) )
    process.hNuMu40puLo = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), systPileupShift = cms.double(-0.7) )

process.hNuEpuHi = process.hNuE.clone( studyMuSelectEff = cms.bool(False), analysisMode = cms.untracked.string('HNUE'), systPileupShift = cms.double(0.84) )
process.hNuEpuLo = process.hNuE.clone( studyMuSelectEff = cms.bool(False), analysisMode = cms.untracked.string('HNUE'), systPileupShift = cms.double(-0.84) )

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

process.hNuMu1QCD = process.hNu.clone( studyMuSelectEff = cms.bool(False), nFakeLeptons = cms.untracked.int32(1) )
process.hNuMu2QCD = process.hNu.clone( studyMuSelectEff = cms.bool(False), nFakeLeptons = cms.untracked.int32(2) )
process.hNuE1QCD  = process.hNuE.clone( studyMuSelectEff = cms.bool(False), nFakeLeptons = cms.untracked.int32(1) )
process.hNuE2QCD  = process.hNuE.clone( studyMuSelectEff = cms.bool(False), nFakeLeptons = cms.untracked.int32(2) )

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
process.hNuTop.studyScaleFactor  = cms.bool(studyTopScale)
process.hNuTop.heepVersion       = cms.int32(heepVersion)
process.hNuTop.electronRho       = cms.InputTag( 'kt6PFJetsForIsolation','rho' )
process.hNuTop.useDirtyElectrons = cms.bool( False ) 
process.hNuTop.useDirtyMuons     = cms.bool( False ) 
process.hNuTop.fakeWeightFlavor  = cms.int32( 0 ) 

#--- Some corrections re-enabled ---#
process.hNuTop.applyMuIDEffcorr      = cms.bool(isMC)
process.hNuTop.applyMuIDEffsign      = cms.int32(0)
process.hNuTop.applyEleEScale        = cms.bool(isMC) 
process.hNuTop.correctEscale         = cms.bool(isMC) 
process.hNuTop.applyEleIDweight      = cms.bool(isMC) 
process.hNuTop.pileupEra             = cms.int32(pileupEra)
process.hNuTop.trigMatchPset.trigEra = cms.int32(dataEra)
process.hNuTop.muIDPset.eraForId     = cms.int32(dataEra)
process.hNuTop.EBidWgt = cms.double( 1.000 ) 
process.hNuTop.EEidWgt = cms.double( 1.000 ) 
#--- Values below zero disable the vertex requirement ---#
process.hNuTop.maxVertexZsepCM     = cms.double(-1)
process.hNuTop.maxVertexJetVZsepCM = cms.double(-1)

process.hNuTop.addSlopeTree        = cms.untracked.bool(addSlopeTrees)

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

process.hNuGoodMuFakeE      = process.hNuTop.clone( useDirtyElectrons = cms.bool(True), useDirtyMuons = cms.bool(False) )
process.hNuFakeMuGoodE      = process.hNuTop.clone( useDirtyElectrons = cms.bool(False), useDirtyMuons = cms.bool(True) )
process.hNuGoodMuFakeEwgtE  = process.hNuTop.clone( useDirtyElectrons = cms.bool(True), useDirtyMuons = cms.bool(False), fakeWeightFlavor = cms.int32(11) )
process.hNuFakeMuGoodEwgtMu = process.hNuTop.clone( useDirtyElectrons = cms.bool(False), useDirtyMuons = cms.bool(True), fakeWeightFlavor = cms.int32(13) )
process.hNuFakeMuEwgtE      = process.hNuTop.clone( useDirtyElectrons = cms.bool(True), useDirtyMuons = cms.bool(True), fakeWeightFlavor = cms.int32(11) )
process.hNuFakeMuEwgtMu     = process.hNuTop.clone( useDirtyElectrons = cms.bool(True), useDirtyMuons = cms.bool(True), fakeWeightFlavor = cms.int32(13) )
process.hNuFakeMuEwgtMuE    = process.hNuTop.clone( useDirtyElectrons = cms.bool(True), useDirtyMuons = cms.bool(True), fakeWeightFlavor = cms.int32(1113) )

#----------------#
#--- LumiList ---#
#----------------#

process.lumilist = cms.EDAnalyzer('LumiList')

if not isMC:
    process.llPath = cms.Path(process.lumilist)

#-------------#
#--- Paths ---#
#-------------#
if runMuonAnalysis:
    if qcdStudy and topStudy:
        process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('pNominal','pMu1QCD','pMu2QCD','pTop') )
    if qcdStudy and not topStudy:
        process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('pNominal','pMu1QCD','pMu2QCD') )
    if topStudy and not qcdStudy:
        process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('pNominal','pTop') )
    if not qcdStudy and not topStudy:
        process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('pNominal') )
else: 
    # if qcdStudy and topStudy:
    #     process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('pMu1QCD','pMu2QCD','pTop') )
    # if qcdStudy and not topStudy:
    #     process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('pMu1QCD','pMu2QCD') )
    # if topStudy and not qcdStudy:
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
                # process.p24mesHi  = cms.Path( process.AnalysisIntroSequence + process.hNuMu24mesHi )
                # process.p24mesLo  = cms.Path( process.AnalysisIntroSequence + process.hNuMu24mesLo )
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
            # process.p40mesHi  = cms.Path( process.AnalysisIntroSequence + process.hNuMu40mesHi )
            # process.p40mesLo  = cms.Path( process.AnalysisIntroSequence + process.hNuMu40mesLo )
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
      process.pEjesHi  = cms.Path(process.AnalysisIntroSequence + process.hNuEjesHi);
      process.pEjesLo  = cms.Path(process.AnalysisIntroSequence + process.hNuEjesLo);
      process.pEidHi   = cms.Path(process.AnalysisIntroSequence + process.hNuEidHi);
      process.pEidLo   = cms.Path(process.AnalysisIntroSequence + process.hNuEidLo);
      process.pEtrigHi = cms.Path(process.AnalysisIntroSequence + process.hNuEtrigHi);
      process.pEtrigLo = cms.Path(process.AnalysisIntroSequence + process.hNuEtrigLo);
      process.pEjerHi  = cms.Path(process.AnalysisIntroSequence + process.hNuEjerHi);
      process.pEjerLo  = cms.Path(process.AnalysisIntroSequence + process.hNuEjerLo);
      process.pEescale = cms.Path(process.AnalysisIntroSequence + process.hNuEescale);
      process.pEpuHi   = cms.Path(process.AnalysisIntroSequence + process.hNuEpuHi);
      process.pEpuLo   = cms.Path(process.AnalysisIntroSequence + process.hNuEpuLo);

   if doTriggerStudy and isData:
       process.pENominal = cms.Path( process.AnalysisIntroSequence + process.TriggerStudyElectronSequence + process.hNuE  )
   else:
       process.pENominal = cms.Path( process.AnalysisIntroSequence + process.hNuE )
          
if qcdStudy:
    process.pMu1QCD = cms.Path( process.AnalysisIntroSequence + process.hNuMu1QCD )
    process.pMu2QCD = cms.Path( process.AnalysisIntroSequence + process.hNuMu2QCD )
    process.pE1QCD  = cms.Path( process.AnalysisIntroSequence + process.hNuE1QCD )
    process.pE2QCD  = cms.Path( process.AnalysisIntroSequence + process.hNuE2QCD )

    process.pGoodMuFakeE      = cms.Path( process.AnalysisIntroSequence + process.hNuGoodMuFakeE )
    process.pFakeMuGoodE      = cms.Path( process.AnalysisIntroSequence + process.hNuFakeMuGoodE )
    process.pGoodMuFakeEwgtE  = cms.Path( process.AnalysisIntroSequence + process.hNuGoodMuFakeEwgtE )
    process.pFakeMuGoodEwgtMu = cms.Path( process.AnalysisIntroSequence + process.hNuFakeMuGoodEwgtMu )
    process.pFakeMuEwgtMu     = cms.Path( process.AnalysisIntroSequence + process.hNuFakeMuEwgtMu )
    process.pFakeMuEwgtE      = cms.Path( process.AnalysisIntroSequence + process.hNuFakeMuEwgtE )
    process.pFakeMuEwgtMuE    = cms.Path( process.AnalysisIntroSequence + process.hNuFakeMuEwgtMuE )

    # if isMC:
    #     if isRun2011A:
    #         process.pQCD24 = cms.Path( process.AnalysisIntroSequence + process.hNuQCDMu24 )
    #     process.pQCD40 = cms.Path( process.AnalysisIntroSequence + process.hNuQCDMu40 ) 
    # else:
    #     process.pQCD = cms.Path( process.AnalysisIntroSequence + process.hNuQCD )
        
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

