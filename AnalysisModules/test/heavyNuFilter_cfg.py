import FWCore.ParameterSet.Config as cms

from operator import isSequenceType
import os

#--- Data/MC switch ---#
isMC=True
isData=not isMC

smode = 0

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

process = cms.Process("PATSKIM");

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)

# source
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('file:/local/cms/user/pastika/heavyNuAnalysis_2012/2C1FBAB2-C1D4-E111-A89A-001E6739815B.root')
)

if isData:
    if isRun2012:
        if isRun2012topoff:
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
            process.GlobalTag.globaltag = cms.string('FT_53_V6_AN1::All')
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
    minJetPt     = cms.double(30.),
    maxMuAbsEta  = cms.double(2.5),
    maxElecAbsEta = cms.double(2.6),
    maxJetAbsEta = cms.double(2.6),
    minMuonJetdR = cms.double(0.3),

    muonTrackRelIsoLimit  = cms.double(0.1), # 10.0),
    maxVertexZsepCM       = cms.double(-1), # disabled
    maxJetVZsepCM         = cms.double(-1), # disabled

    isPFJets = cms.bool(True),

    heepVersion = cms.untracked.int32(-1),

    mode = cms.int32(smode), #0 for Z+Jets, 1 for TTBar and other
    electronRho  = cms.InputTag( 'kt6PFJetsForIsolation','rho' )
    )

process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') )
process.out.outputCommands = cms.untracked.vstring("drop *", "keep *_*_*_RECO", "keep *_*_*_LHE", "keep *_*_*_SIM", "keep *_*_*_HLT", "keep *_*_*_PATSKIM")
#                                                   "drop patTaus_*_*_PATSKIM", "drop recoPFCandidates_*_*_PATSKIM", "drop recoIsoDepositedmValueMap_*_*_PATSKIM")
process.p = cms.Path( process.AnalysisIntroSequence * process.hNu )
process.end = cms.EndPath(process.out)
