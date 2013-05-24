import FWCore.ParameterSet.Config as cms

from operator import isSequenceType
import os, random

#--- Data/MC switch ---#
isMC=False
isData=not isMC

#--- Specify CMSSW release (53, 52, 44, 42, 37) ---#
cmsswRelease = 53

#--- Signal MC flags ---#
isMCsignal=False

#--- Data Run era flags ---#
# options 2012AB, 2012Ar, 2012Cr, 2012Cp, 2012Ce, 2012D
# options 2012ABCre, 2012Dre
runEra = "2012ABCre"

#--- Flags for data taking era, which are set automatically ---#
#--- dataEra has been deprecated ---#
# dataEra   = 20122 
pileupEra = 20122

#--- Filtering for multi-skims ---#
lljjAnalysisSkim = False
llPrescaleSkim   = False
llHighMassSkim   = False
topSkim          = False

#--- Flags for nominal studies ---#
runMuonAnalysis     = False
runElectronAnalysis = True
systematics    = False
tagandprobe    = False
doTriggerStudy = False
addSlopeTrees  = True

#--- HEEP ID for electrons ---#
#--- Recognized values: 41 or 40 (2012), 31 or 32 (2011) ---#
heepVersion = 41

#--- Flags for Top studies ---#
topStudy  = True

#--- Flags for QCD studies ---#
qcdStudy  = False
doDijet   = False

if isData:
    systematics      = False
    qcdStudy         = False
    lljjAnalysisSkim = False
    if runMuonAnalysis:
        runElectronAnalysis = False
        topStudy = True
        topSkim  = False
    if runElectronAnalysis:
        runMuonAnalysis = False
        topSkim  = False
        topStudy = False

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
                            fileNames=cms.untracked.vstring('file:/local/cms/user/pastika/heavyNuAnalysis_2012/skims/2012A_Jan22_2013_sample.root')
                            #file:/local/cms/user/pastika/heavyNuAnalysis_2012/skims/sherpadyfile.root')
                            #file:/hdfs/cms/skim/mu/hNu_2012/skim50_35/jul13rereco_run2012AB/muTopMultiSkim_plep50_35_pemu50_35_2012ABjul13_feb11_104.root')
                            #/hdfs/cms/phedex/store/mc/Summer12_DR53X/WRToNuLeptonToLLJJ_MW-2900_MNu-1450_TuneZ2star_8TeV-pythia6-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/30672247-B1EC-E111-A906-00215E222220.root')
                            #/store/mc/Summer12_DR53X/WRToNuLeptonToLLJJ_MW-2900_MNu-1450_TuneZ2star_8TeV-pythia6-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/30672247-B1EC-E111-A906-00215E222220.root')
                            #/store/mc/Summer12_DR53X/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/AODSIM/PU_S10_START53_V7C-v1/00000/56B5A94D-FB27-E211-901F-AC162DA8C2B0.root')
                            #file:/local/cms/user/pastika/heavyNuAnalysis_2012/skims/DY0JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_START53_V7A/FZJ_000.root')
                            #file:/local/cms/user/pastika/heavyNuAnalysis_2012/skims/TTBar_Partial/heavynu_candevents_425_1_uct.root')
                            #file:/local/cms/user/pastika/heavyNuAnalysis_2012/skims/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_START53_V7A_2/heavynu_candevents_136_1_Hd1.root')
                            #file:/hdfs/cms/skim/elec/hNu_2012/photon/eejj_skim_aug8_rereco_run2012a/eejj_skim_aug8_rereco_run2012a_266.root')
                            #file:/hdfs/cms/skim/mu/hNu_2012/jul13_2012A_mu35e35/jul13_2012A_mu35e35_015.root')
                            #file:/local/cms/user/pastika/heavyNuAnalysis_2012/2C1FBAB2-C1D4-E111-A89A-001E6739815B.root
                            #file:/home/ugrad/pastika/cms/HeavyNu/CMSSW_5_3_3_patch1/src/HeavyNu/AnalysisModules/heavynu_candevents.root')
                            #file:/hdfs/cms/skim/mu/hNu_2012/jul13_2012A_mu35e35/jul13_2012A_mu35e35_015.root
)

if isData:
    if runEra == "2012AB" or runEra == "2012Ar" or runEra == "2012Cr" or runEra == "2012Ce":
        from HeavyNu.AnalysisModules.goodLumiList_2012abc_rereco_cfi import lumisToProcess
    elif runEra == "2012Cp" or runEra == "2012D":
        from HeavyNu.AnalysisModules.goodLumiList_2012cd_prompt_cfi import lumisToProcess
    elif runEra == "2012ABCre" or runEra == "2012Dre":
        from HeavyNu.AnalysisModules.goodLumiList_2012abcd_rereco_22jan2013_cfi import lumisToProcess
                            
    #process.source.lumisToProcess = lumisToProcess

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


## Load additional processes
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

## Global Tags:
if (isMC):
    print "=================> MC flag is SET <===================="
    if cmsswRelease == 52:
        process.GlobalTag.globaltag=cms.string('START52_V11::All')
    elif cmsswRelease == 53:
        process.GlobalTag.globaltag=cms.string('START53_V15::All')
    else:
        print "INVALID CMSSW release id %(rid)i"%{"rid":cmsswRelease}
else:
    print "===============> Running on DATA <===================="
    if cmsswRelease == 52:
        process.GlobalTag.globaltag = cms.string('GR_R_52_V9::All')
    elif cmsswRelease == 53:
        if runEra == "2012AB":
            process.GlobalTag.globaltag = cms.string('FT_53_V6_AN3::All')
        elif runEra == "2012Ar":
            process.GlobalTag.globaltag = cms.string('FT_53_V6C_AN3::All') 
        elif runEra == "2012Cr":
            process.GlobalTag.globaltag = cms.string('FT_53_V10_AN3::All')
        elif runEra == "2012Cp":
            process.GlobalTag.globaltag = cms.string('GR_P_V41_AN3::All')
        elif runEra == "2012Ce":
            process.GlobalTag.globaltag = cms.string('FT_P_V42C_AN3::All')
        elif runEra == "2012D":
            process.GlobalTag.globaltag = cms.string('GR_P_V42_AN3::All')
        elif runEra == "2012ABCre":
            process.GlobalTag.globaltag = cms.string('FT_R_53_V18::All')
        elif runEra == "2012Dre":
            process.GlobalTag.globaltag = cms.string('FT_R_53_V21::All')
    else:
        print "INVALID CMSSW release id %(rid)i"%{"rid":cmsswRelease}

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')

################################################################################################
###   P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D   ###
################################################################################################

## pat sequences to be loaded:
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.pfTools import *

def usePF2PAT_WREdition(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix="", jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute']), pvCollection=cms.InputTag('offlinePrimaryVertices'), typeIMetCorrections=False, outputModules=['out']):
    # PLEASE DO NOT CLOBBER THIS FUNCTION WITH CODE SPECIFIC TO A GIVEN PHYSICS OBJECT.
    # CREATE ADDITIONAL FUNCTIONS IF NEEDED.

    """Switch PAT to use PF2PAT instead of AOD sources. if 'runPF2PAT' is true, we'll also add PF2PAT in front of the PAT sequence"""

    # -------- CORE ---------------
    if runPF2PAT:
        process.load("CommonTools.ParticleFlow.PF2PAT_cff")
        process.patPF2PATSequence = cms.Sequence( process.PF2PAT + process.patDefaultSequence)
    else:
        process.patPF2PATSequence = cms.Sequence( process.patDefaultSequence )

    if not postfix == "":
        from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet
        cloneProcessingSnippet(process, process.patPF2PATSequence, postfix)

    removeCleaning(process, postfix=postfix, outputModules=outputModules)

    # Photons
    print "Temporarily switching off photons completely"

    removeSpecificPATObjects(process,names=['Photons', 'Taus', "Muons", "Electrons"],outputModules=outputModules,postfix=postfix)
    removeIfInSequence(process,"patPhotonIsolation","patDefaultSequence",postfix)

    # Jets
    if runOnMC :
        switchToPFJets( process, cms.InputTag('pfNoTau'+postfix), jetAlgo, postfix=postfix,
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
        switchToPFJets( process, cms.InputTag('pfNoTau'+postfix), jetAlgo, postfix=postfix,
                        jetCorrections=jetCorrections, type1=typeIMetCorrections, outputModules=outputModules )
                        

    # Taus
    if ishpsPFTau:
        adaptPFTaus( process, tauType='hpsPFTau', postfix=postfix )

    # MET
    switchToPFMET(process, cms.InputTag('pfMET'+postfix), type1=typeIMetCorrections, postfix=postfix)
    if not runOnMC :
        if hasattr(process,'patPFMet'+postfix):
            getattr(process,'patPFMet'+postfix).addGenMET = cms.bool(False)

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
elif cmsswRelease == 53:
	if isMC:
		process.GlobalTag.toGet = cms.VPSet(
		  cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
		       tag = cms.string("TrackProbabilityCalibration_2D_Data53X_v2"),
		       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
		  cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
		       tag = cms.string("TrackProbabilityCalibration_3D_Data53X_v2"),
		       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
		)
	else:
		process.GlobalTag.toGet = cms.VPSet(
		  cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
		       tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
		       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
		  cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
		     tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
		       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
		)


## --------------------- ##
## Define the basic path ##
## --------------------- ##

#--- Trigger filters for new skims ---#
process.pickTriggered = cms.EDFilter('HeavyNuSkimFilter',
    triggerResults = cms.InputTag('TriggerResults','','HNUSKIMS'), 
    filterOnPaths  = cms.vstring('plep35_35_prescale')  
)
process.twoLeptonSkimFilter         = process.pickTriggered.clone( filterOnPaths = cms.vstring('plep50_35') )
process.twoLeptonPrescaleSkimFilter = process.pickTriggered.clone( filterOnPaths = cms.vstring('plep35_35_prescale') )
process.twoLeptonHighMassSkimFilter = process.pickTriggered.clone( filterOnPaths = cms.vstring('plep35_35_mLL170') )
process.topSkimFilter               = process.pickTriggered.clone( filterOnPaths = cms.vstring('pemu50_35') )
process.topAndTwoLeptonSkimFilter   = process.pickTriggered.clone( filterOnPaths = cms.vstring('plep50_35','pemu50_35') )

process.boolTrue = cms.EDFilter( 'HLTBool',
    result = cms.bool( True )
)
process.skimFilterSequence = cms.Sequence(process.boolTrue)
if lljjAnalysisSkim:
    if topSkim:
        process.skimFilterSequence += process.topAndTwoLeptonSkimFilter
    else:
        process.skimFilterSequence += process.twoLeptonSkimFilter
if llPrescaleSkim:
    process.skimFilterSequence += process.twoLeptonPrescaleSkimFilter
if llHighMassSkim:
    process.skimFilterSequence += process.twoLeptonHighMassSkimFilter

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

#--- ECAL laser noise filter ---#
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')

process.eventFilters = cms.Sequence( process.skimFilterSequence * process.scrapingFilter * process.primaryVertexFilter * process.HBHENoiseFilter * process.ecalLaserCorrFilter) 

#--- HCAL laser events filter ---#
if isData:
    process.load("EventFilter.HcalRawToDigi.hcallasereventfilter2012_cff")
    process.eventFilters += process.hcallLaserEvent2012Filter

#--- Pre-filter ---#
process.prefilter = cms.EDFilter("HeavyNuPreFilter",
                                 muonTag      = cms.InputTag( 'selectedPatMuons' ),
                                 electronTag  = cms.InputTag( 'selectedPatElectrons' ),
                                 electronRho  = cms.InputTag( 'kt6PFJetsForIsolation','rho' ),
                                 ptL1         = cms.untracked.double(40.0),
                                 ptL2         = cms.untracked.double(40.0))

if isMC:
   # Gen Level Energy balance filter to fix Pythia6 lhe interface bug
   process.load("HeavyNu.AnalysisModules.hnuTotalKinematicsFilter_cfi")
   process.AnalysisIntroSequence = cms.Sequence(
#       process.hnuTotalKinematicsFilter * process.eventFilters * process.patDefaultSequence * process.patTrackSequence * process.kt6PFJetsForIsolation
       process.eventFilters * process.patDefaultSequence * process.patTrackSequence * process.kt6PFJetsForIsolation * process.prefilter
   )
else:
   process.AnalysisIntroSequence = cms.Sequence(
       process.eventFilters * process.patDefaultSequence * process.patTrackSequence * process.kt6PFJetsForIsolation * process.prefilter
   )
if isPFJets:
    process.AnalysisIntroSequence *= process.modifiedPF2PATSequence

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
    removeSpecificPATObjects(process, names = ['Jets','Taus'], outputModules = [])

# Special change for saving good products in data
if isData:
    process.out.outputCommands = cms.untracked.vstring("keep *")

#--- Trigger matching ---#
process.load("HeavyNu.AnalysisModules.hnutrigmatch_cfi")
from PhysicsTools.PatAlgos.tools.trigTools import *
if isData:
    if (runMuonAnalysis or topStudy) and runElectronAnalysis:
        triggerMatchersList = [ 'muonTriggerMatchHLTMuons', 'electronTriggerMatchHLTElectrons' ]
    elif runMuonAnalysis or topStudy:
        triggerMatchersList = [ 'muonTriggerMatchHLTMuons', ]
    elif runElectronAnalysis:
        triggerMatchersList = [ 'electronTriggerMatchHLTElectrons', ]
	switchOnTriggerMatching( process, triggerMatchers = triggerMatchersList, outputModule = '' )
    
    if runMuonAnalysis or topStudy:
        switchOnTriggerMatchEmbedding( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ], outputModule = '' )
        removeCleaningFromTriggerMatching( process, outputModule = '' )
        # if runEra == "2012AB" or "2012Cr" or "2012Cp" or "2012D":
        process.muonTriggerMatchHLTMuons.matchedCuts = cms.string( 'path( "HLT_Mu40_eta2p1_v*",1,0 )' )
    if runElectronAnalysis:
        switchOnTriggerMatchEmbedding( process, triggerMatchers = [ 'electronTriggerMatchHLTElectrons' ], outputModule = '' )
        removeCleaningFromTriggerMatching( process, outputModule = '' )
        # if runEra == "2012AB" or "2012Cr" or "2012Cp" or "2012D":
        process.electronTriggerMatchHLTElectrons.matchedCuts = cms.string( 'path( "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*",1,0 )' )


#--- Electrons trigger analysis ---#
#process.dumpEvContent = cms.EDAnalyzer("EventContentAnalyzer")
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
#if isMC and isMCsignal:
#    process.hNu.isSignal = cms.bool(True)
#    process.load("HeavyNu.AnalysisModules.heavyNuGenFilter_cfi")
#    process.hNuGenFilter.keepIds = cms.vint32(2,)

process.hNu.minMu2pt         = cms.double(40.)
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
process.hNu.randseed     = cms.untracked.uint32(random.randint(0, 1 << 32 -1))

#--- Muon ID corrections are taken from June 22, 2012 studies---#
process.hNu.applyMuIDEffcorr = cms.bool(isMC)

if isData:
    if runMuonAnalysis:
        process.hNu.muonTag     = cms.InputTag( 'selectedPatMuonsTriggerMatch' )
    else:
        process.hNu.muonTag     = cms.InputTag( 'selectedPatMuons' )
    process.hNu.trigMatchPset.trigEventTag    = cms.InputTag("patTriggerEvent")
    process.hNu.trigMatchPset.electronFilters = cms.vstring('hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter')
    process.hNu.muIDPset.eraForId             = cms.int32( 1 )
    process.hNu.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_v1','HLT_Mu40_v2','HLT_Mu40_v3','HLT_Mu40_v5','HLT_Mu40_eta2p1_v1','HLT_Mu40_eta2p1_v2','HLT_Mu40_eta2p1_v3','HLT_Mu40_eta2p1_v4','HLT_Mu40_eta2p1_v5','HLT_Mu40_eta2p1_v6','HLT_Mu40_eta2p1_v7','HLT_Mu40_eta2p1_v8','HLT_Mu40_eta2p1_v9','HLT_Mu40_eta2p1_v10','HLT_Mu40_eta2p1_v11' ) 
    process.hNu.trigMatchPset.electronTriggers = cms.vstring( 'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3','HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v4','HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v5','HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v6','HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7' )
    process.hNu.trigMatchPset.triggerPt    = cms.double( 40. )
    process.hNu.trigMatchPset.muTriggerEta = cms.double( 2.1 )
else:
    # turn on MC trigger simulation
    process.hNu.trigMatchPset.randomSeed = cms.int32(os.getpid())
    # process.hNu.trigMatchPset.trigEra    = cms.int32( dataEra )
    # Parameters for muon ID corrections to MC
    process.hNu.muIDPset.eraForId        = cms.int32( 1 )
    process.hNu.trigMatchPset.triggerPt    = cms.double( 40. )
    process.hNu.trigMatchPset.muTriggerEta = cms.double( 2.1 )

if isPFJets:
    process.hNu.jetTag = cms.InputTag( 'selectedPatJetsPFlow' )

# heavyNuGenFilter setup
process.load("HeavyNu.AnalysisModules.heavyNuGenFilter_cfi")
#--- Necessary to sort out trigger complications ---#
process.hNuMu40       = process.hNu.clone() 
process.hNuMu40eta2p1 = process.hNu.clone()

process.hNuE               = process.hNu.clone(analysisMode = cms.untracked.string('HNUE'))
process.hNuE.correctEscale = cms.bool(isMC)

process.hNuEMu               = process.hNu.clone(analysisMode = cms.untracked.string('TOP'))
process.hNuEMu.correctEscale = cms.bool(isMC)

process.hTauX               = process.hNu.clone(analysisMode = cms.untracked.string('TAUX'))

process.hNuETrig = process.hNuE.clone()

# signal gen filter
process.load("HeavyNu.AnalysisModules.heavyNuGenFilter_cfi")
if isMC and isMCsignal:
	process.hNuGenFilterElec = process.hNuGenFilter.clone(keepIds = cms.vint32(1,))
	process.hNuGenFilterMu   = process.hNuGenFilter.clone(keepIds = cms.vint32(2,))
	process.hNuGenFilterTau  = process.hNuGenFilter.clone(keepIds = cms.vint32(3,))

#--- Electron ID corrections are taken from June 22, 2012 studies---#
process.hNu.applyMuIDEffcorr = cms.bool(isMC)

process.hNuMu40.trigMatchPset.triggerPt    = cms.double( 40. )
process.hNuMu40.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_v1','HLT_Mu40_v2','HLT_Mu40_v3','HLT_Mu40_v5' ) 
process.hNuMu40.trigMatchPset.randomSeed   = cms.int32( os.getpid() )
process.hNuMu40eta2p1.trigMatchPset.triggerPt    = cms.double( 40. )
process.hNuMu40eta2p1.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_eta2p1_v1','HLT_Mu40_eta2p1_v2','HLT_Mu40_eta2p1_v3','HLT_Mu40_eta2p1_v4','HLT_Mu40_eta2p1_v5','HLT_Mu40_eta2p1_v6','HLT_Mu40_eta2p1_v7','HLT_Mu40_eta2p1_v8','HLT_Mu40_eta2p1_v9' ) 
process.hNuMu40eta2p1.trigMatchPset.randomSeed   = cms.int32( os.getpid() )

process.hNuMu40jesHi = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyJECUsign = cms.int32(1) )
process.hNuMu40jesLo = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyJECUsign = cms.int32(-1) )

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


process.hNuMu40mer = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), checkMERUnc = cms.bool(True) )

process.hNuMu40midHi = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyMuIDEffsign = cms.int32(1) )
process.hNuMu40midLo = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyMuIDEffsign = cms.int32(-1) )
process.hNuMu40midHi.applyMuIDEffcorr = cms.bool( isMC )
process.hNuMu40midLo.applyMuIDEffcorr = cms.bool( isMC )

process.hNuMu40trigHi = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyTrigEffsign  = cms.int32(1) )
process.hNuMu40trigLo = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), applyTrigEffsign  = cms.int32(-1) )

# Pileup uncertainty: +/- 8% on number of interactions leads to 0.4 in 2011A, 0.7 in 2011B
#                     +/- 5% on number of interactions (16.85, 12/06/09) leads to 0.84 in 2012AB
if runEra == "2012AB" or "2012Cr" or "2012Cp" or "2012D":
    process.hNuMu40puHi = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), systPileupShift = cms.int32(1) )
    process.hNuMu40puLo = process.hNuMu40.clone( studyMuSelectEff = cms.bool(False), systPileupShift = cms.int32(-1) )

process.hNuEpuHi = process.hNuE.clone( studyMuSelectEff = cms.bool(False), analysisMode = cms.untracked.string('HNUE'), systPileupShift = cms.int32(1) )
process.hNuEpuLo = process.hNuE.clone( studyMuSelectEff = cms.bool(False), analysisMode = cms.untracked.string('HNUE'), systPileupShift = cms.int32(-1) )

## ============ ##
## QCD Analysis ##
## ============ ##
process.load("HeavyNu.AnalysisModules.heavynuqcd_cfi")

process.hNuQCD.isPFJets = cms.bool(isPFJets)
process.hNuQCD.getSurvivalRate = cms.bool( doDijet )    
#--- Pileup corrections: needed even for data! ---#
process.hNuQCD.pileupEra = cms.int32(pileupEra)

#--- Muon QCD estimates and closure tests ---#
process.hNuMu1QCD = process.hNu.clone( studyMuSelectEff = cms.bool(False), nFakeLeptons = cms.untracked.int32(1) )
process.hNuMu2QCD = process.hNu.clone( studyMuSelectEff = cms.bool(False), nFakeLeptons = cms.untracked.int32(2) )
#--- Fakes start with "loose" muons ---#
process.hNuMu1QCD.muIDPset.eraForId = cms.int32( -1 )
process.hNuMu2QCD.muIDPset.eraForId = cms.int32( -1 )

process.hNuMu1QCDMiJ = process.hNuMu1QCD.clone( qcdMuonsInJets = cms.untracked.bool(True), qcdMuJetOverlap = cms.double(0.5) ) 
process.hNuMu2QCDMiJ = process.hNuMu2QCD.clone( qcdMuonsInJets = cms.untracked.bool(True), qcdMuJetOverlap = cms.double(0.5) ) 
process.hNuMu1QCDMiJ0p3 = process.hNuMu1QCD.clone( qcdMuonsInJets = cms.untracked.bool(True), qcdMuJetOverlap = cms.double(0.3) ) 
process.hNuMu2QCDMiJ0p3 = process.hNuMu2QCD.clone( qcdMuonsInJets = cms.untracked.bool(True), qcdMuJetOverlap = cms.double(0.3) ) 

#--- Electron QCD estimates and closure tests ---#
process.hNuE1QCD = process.hNuE.clone( studyMuSelectEff = cms.bool(False), nFakeLeptons = cms.untracked.int32(1) )
process.hNuE2QCD = process.hNuE.clone( studyMuSelectEff = cms.bool(False), nFakeLeptons = cms.untracked.int32(2) )

if isData:
    process.hNuQCD.muonTag                     = cms.InputTag( 'selectedPatMuonsTriggerMatch' )
    process.hNuQCD.trigMatchPset.trigEventTag  = cms.InputTag( 'patTriggerEvent' )
    process.hNuQCD.trigMatchPset.electronFilters  = cms.vstring('') 
    process.hNuQCD.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_v1','HLT_Mu40_v2','HLT_Mu40_v3','HLT_Mu40_v5','HLT_Mu40_eta2p1_v1','HLT_Mu40_eta2p1_v2','HLT_Mu40_eta2p1_v3','HLT_Mu40_eta2p1_v4','HLT_Mu40_eta2p1_v5','HLT_Mu40_eta2p1_v6','HLT_Mu40_eta2p1_v7','HLT_Mu40_eta2p1_v8','HLT_Mu40_eta2p1_v9','HLT_Mu40_eta2p1_v10','HLT_Mu40_eta2p1_v11','HLT_Mu40_eta2p1_v12' ) 
    process.hNuQCD.trigMatchPset.triggerPt    = cms.double( 40. )
    process.hNuQCD.trigMatchPset.electronTriggers = cms.vstring( '' )
else:
    # turn on MC trigger simulation
    process.hNuQCD.trigMatchPset.randomSeed = cms.int32( os.getpid() )
    # process.hNuQCD.trigMatchPset.trigEra    = cms.int32( dataEra )
    # Parameters for muon ID corrections
    process.hNuQCD.muIDPset.eraForId        = cms.int32( 1 )

if isPFJets:
    process.hNuQCD.jetTag  = cms.InputTag( 'selectedPatJetsPFlow' )
    if doDijet:
        process.hNuQCD.metTag  = cms.InputTag( 'pfType1CorrectedMet' )

process.hNuQCDMu40eta2p1 = process.hNuQCD.clone() 
process.hNuQCDMu40eta2p1.trigMatchPset.triggerPt  = cms.double( 40. )
process.hNuQCDMu40eta2p1.trigMatchPset.randomSeed = cms.int32( os.getpid() )

process.hNuQCDMu40 = process.hNuQCD.clone() 
process.hNuQCDMu40.trigMatchPset.triggerPt  = cms.double( 40. )
process.hNuQCDMu40.trigMatchPset.randomSeed = cms.int32( os.getpid() )
process.hNuQCDMu40.trigMatchPset.muonTriggers = cms.vstring( 'HLT_Mu40_v1','HLT_Mu40_v2','HLT_Mu40_v3','HLT_Mu40_v4','HLT_Mu40_v5','HLT_Mu40_v6','HLT_Mu40_v7','HLT_Mu40_v8','HLT_Mu40_v9','HLT_Mu40_v10','HLT_Mu40_v11','HLT_Mu40_v12','HLT_Mu40_v13','HLT_Mu40_v14','HLT_Mu40_v15','HLT_Mu40_v16','HLT_Mu40_v17','HLT_Mu40_v18','HLT_Mu40_v19','HLT_Mu40_v20' ) 

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
            process.AnalysisIntroSequenceMu = cms.Sequence( process.hNuGenFilterMu * process.AnalysisIntroSequence )
        else:
            process.AnalysisIntroSequenceMu = process.AnalysisIntroSequence

        process.p40 = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40 ) 

        if systematics:

            process.p40jesHi  = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40jesHi )
            process.p40jesLo  = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40jesLo )
            process.p40jerHi  = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40jerHi )
            process.p40jerLo  = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40jerLo )
            # process.p40mesHi  = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40mesHi )
            # process.p40mesLo  = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40mesLo )
            process.p40mer    = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40mer )
            process.p40midHi  = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40midHi )
            process.p40midLo  = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40midLo )
            process.p40trigHi = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40trigHi )
            process.p40trigLo = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40trigLo )
            process.p40puHi   = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40puHi )
            process.p40puLo   = cms.Path( process.AnalysisIntroSequenceMu + process.hNuMu40puLo )
    else:
        process.pNominal = cms.Path( process.AnalysisIntroSequence + process.hNu )
        if doTriggerStudy:
            # if runEra == "2012AB" or "2012C" or "2012D":
            process.p40eta2p1 = cms.Path( process.AnalysisIntroSequence + process.hNuMu40eta2p1 )

if runElectronAnalysis:
    if isMC:
        if isMCsignal:
            process.AnalysisIntroSequenceElec = cms.Sequence( process.hNuGenFilterElec * process.AnalysisIntroSequence )
        else:
            process.AnalysisIntroSequenceElec = process.AnalysisIntroSequence

        process.pE = cms.Path(process.AnalysisIntroSequenceElec + process.hNuE)
        if systematics:
            process.pEjesHi  = cms.Path(process.AnalysisIntroSequenceElec + process.hNuEjesHi);
            process.pEjesLo  = cms.Path(process.AnalysisIntroSequenceElec + process.hNuEjesLo);
            process.pEidHi   = cms.Path(process.AnalysisIntroSequenceElec + process.hNuEidHi);
            process.pEidLo   = cms.Path(process.AnalysisIntroSequenceElec + process.hNuEidLo);
            process.pEtrigHi = cms.Path(process.AnalysisIntroSequenceElec + process.hNuEtrigHi);
            process.pEtrigLo = cms.Path(process.AnalysisIntroSequenceElec + process.hNuEtrigLo);
            process.pEjerHi  = cms.Path(process.AnalysisIntroSequenceElec + process.hNuEjerHi);
            process.pEjerLo  = cms.Path(process.AnalysisIntroSequenceElec + process.hNuEjerLo);
            process.pEescale = cms.Path(process.AnalysisIntroSequenceElec + process.hNuEescale);
            process.pEpuHi   = cms.Path(process.AnalysisIntroSequenceElec + process.hNuEpuHi);
            process.pEpuLo   = cms.Path(process.AnalysisIntroSequenceElec + process.hNuEpuLo);

    if isData:
        process.pE = cms.Path( process.AnalysisIntroSequence + process.hNuE )
        if doTriggerStudy:
            process.pE = cms.Path( process.AnalysisIntroSequence + process.TriggerStudyElectronSequence + process.hNuETrig  )
          
if qcdStudy:
    if runMuonAnalysis:
        process.pMu1QCD = cms.Path( process.AnalysisIntroSequence + process.hNuMu1QCD )
        process.pMu2QCD = cms.Path( process.AnalysisIntroSequence + process.hNuMu2QCD )
        process.pMu1QCDMiJ = cms.Path( process.AnalysisIntroSequence + process.hNuMu1QCDMiJ )
        process.pMu2QCDMiJ = cms.Path( process.AnalysisIntroSequence + process.hNuMu2QCDMiJ )
        process.pMu1QCDMiJ0p3 = cms.Path( process.AnalysisIntroSequence + process.hNuMu1QCDMiJ0p3 )
        process.pMu2QCDMiJ0p3 = cms.Path( process.AnalysisIntroSequence + process.hNuMu2QCDMiJ0p3 )
    if runElectronAnalysis:
        process.pE1QCD  = cms.Path( process.AnalysisIntroSequence + process.hNuE1QCD )
        process.pE2QCD  = cms.Path( process.AnalysisIntroSequence + process.hNuE2QCD )

    if doDijet:
        process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

        process.pDijet = cms.Path( process.AnalysisIntroSequence + process.producePFMETCorrections + process.hNuQCD )
        if isData:
            process.pDijetMu40       = cms.Path( process.AnalysisIntroSequence + process.producePFMETCorrections + process.hNuQCDMu40 )
            process.pDijetMu40eta2p1 = cms.Path( process.AnalysisIntroSequence + process.producePFMETCorrections + process.hNuQCDMu40eta2p1 )

    # process.pGoodMuFakeE      = cms.Path( process.AnalysisIntroSequence + process.hNuGoodMuFakeE )
    # process.pFakeMuGoodE      = cms.Path( process.AnalysisIntroSequence + process.hNuFakeMuGoodE )
    # process.pGoodMuFakeEwgtE  = cms.Path( process.AnalysisIntroSequence + process.hNuGoodMuFakeEwgtE )
    # process.pFakeMuGoodEwgtMu = cms.Path( process.AnalysisIntroSequence + process.hNuFakeMuGoodEwgtMu )
    # process.pFakeMuEwgtMu     = cms.Path( process.AnalysisIntroSequence + process.hNuFakeMuEwgtMu )
    # process.pFakeMuEwgtE      = cms.Path( process.AnalysisIntroSequence + process.hNuFakeMuEwgtE )
    # process.pFakeMuEwgtMuE    = cms.Path( process.AnalysisIntroSequence + process.hNuFakeMuEwgtMuE )
        
if topStudy:
    if isMC:
        process.pEMU = cms.Path( process.AnalysisIntroSequence + process.hNuEMu )
    else:
        process.pEMu = cms.Path( process.AnalysisIntroSequence + process.hNuEMu )

if isMCsignal and isMC:
	process.pTauX = cms.Path( process.hNuGenFilterTau * process.AnalysisIntroSequence + process.hTauX )
	

