import FWCore.ParameterSet.Config as cms

from operator import isSequenceType
import os

#--- Data/MC switch ---#
isMC=False
isData=not isMC

#--- Specify CMSSW release (53 is all that is considered for skimming) ---#
cmsswRelease = 53

#--- Signal MC flags ---#
isMCsignal=False

#--- Data Run era flags ---#
# options 2012AB, 2012Ar, 2012Cr, 2012Cp, 2012D
runEra = "2012Cp"

#--- Skim flags ---#
runMuonSkim     = True
runElectronSkim = False
runTopSkim      = True

#--- Prescale for basic path ---#
thePrescaleFactor = 20

if runElectronSkim:
    runTopSkim = False

#--- Tau Flag ---#
ishpsPFTau = False

#--- Should always be True ---#
isPileupMC = True
isPFJets   = True

process = cms.Process("HNUSKIMS");

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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


## Load additional processes
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

## Global Tags:
if (isMC):
    print "=================> MC flag is SET <===================="
    if cmsswRelease == 53:
        process.GlobalTag.globaltag=cms.string('START53_V15::All')
    else:
        print "INVALID CMSSW release id %(rid)i"%{"rid":cmsswRelease}
else:
    print "===============> Running on DATA <===================="
    if cmsswRelease == 53:
        if runEra == "2012AB":
            process.GlobalTag.globaltag = cms.string('FT_53_V6_AN3::All')
        elif runEra == "2012Ar":
            process.GlobalTag.globaltag = cms.string('FT_53_V6C_AN3::All') 
        elif runEra == "2012Cr":
            process.GlobalTag.globaltag = cms.string('FT_53_V10_AN3::All')
        elif runEra == "2012Cp":
            process.GlobalTag.globaltag = cms.string('GR_P_V41_AN3::All')
        elif runEra == "2012D":
            process.GlobalTag.globaltag = cms.string('GR_P_V42_AN3::All')
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
    maxSize = cms.untracked.int32(3000000),
    # save only events passing the full path
    # SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('plep35_35_prescale','plep50_35','plep35_35_mLL170','pemu50_35')),
    SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('plep35_35_prescale','plep50_35','pemu50_35')),
    outputCommands = cms.untracked.vstring("keep *","drop *_*_*_HNUSKIMS","keep edmTriggerResults_TriggerResults__HNUSKIMS")
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
    cut = cms.string('pt > 20')
)
from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import patGenericParticles
process.allPatTracks = patGenericParticles.clone(
    src = cms.InputTag("patAODTrackCands")
)
from PhysicsTools.PatAlgos.selectionLayer1.trackSelector_cfi import *
process.patTracksPt30 = selectedPatTracks.clone(
    cut = 'pt > 30.'
)
process.patTrackSequence = cms.Sequence( 
        process.patAODTrackCandsUnfiltered *
        process.patAODTrackCands *
        process.allPatTracks *
        process.patTracksPt30
)

## ---- ##
## Btag ##
## ---- ##

## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagJetProbabilityCalibration?redirectedfrom=CMS.SWGuideBTagJetProbabilityCalibration#Calibration_in_52x_and_53x_Data

if cmsswRelease == 53:
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


## ----------------------------------- ##
## Define the basic path for each skim ##
## ----------------------------------- ##
process.patSkimCandidateSummary = process.patCandidateSummary.clone() 
process.patSkimCandidateSummary.candidates = cms.VInputTag( cms.InputTag("patMuons"),cms.InputTag("patElectrons") )
process.patSkimCandidates = cms.Sequence( process.makePatMuons + process.makePatElectrons + process.patSkimCandidateSummary )
process.patDefaultSequence = cms.Sequence( process.patSkimCandidates )

process.boolTrue = cms.EDFilter( 'HLTBool',
    result = cms.bool( True )
)
process.boolFalse = cms.EDFilter( 'HLTBool',
    result = cms.bool( False )
)
process.preSkim = cms.Sequence( process.boolTrue * process.patDefaultSequence * process.patTrackSequence * process.kt6PFJetsForIsolation ) 

#--- Muons and top-like ---#
process.muPreSkim = cms.EDFilter("MuFilter",
    muonTag  = cms.InputTag("patMuons"),
    trackTag = cms.InputTag("patTracksPt30"),
    ebTag    = cms.InputTag("correctedHybridSuperClusters"),
    eeTag    = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),

    minMuonPt = cms.double( 35.0 ),
    minSCEt   = cms.double( 35.0 ), 
    overlap   = cms.double( 0.05 ),
    trackPrescale = cms.int32( -1 ),
    trackOnly = cms.bool(False),
    trackMassLow = cms.double( 40.0),
    trackMassHigh = cms.double( 150.0)
)

#--- Electrons ---#
process.electronsAbove35 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("patElectrons"),
    cut = cms.string('(ecalEnergy*sin(superClusterPosition.theta)) >' + str(35) )
)
process.twoElectronsAbove35 = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("ElectronsAbove35"),
    minNumber = cms.uint32(2)
)
process.elecPreSkim = cms.Sequence( process.electronsAbove35 + process.twoElectronsAbove35 )

process.plep35_35_prescale = cms.Path( process.preSkim )
process.plep50_35          = cms.Path( process.preSkim )
process.plep35_35_mLL170   = cms.Path( process.preSkim )

process.pemu50_35          = cms.Path( process.preSkim )

if runMuonSkim:
    process.plep35_35_prescale += process.muPreSkim
    process.plep50_35          += process.muPreSkim
    process.plep35_35_mLL170   += process.muPreSkim
elif runElectronSkim:
    process.plep35_35_prescale += process.elecPreSkim
    process.plep50_35          += process.elecPreSkim
    process.plep35_35_mLL170   += process.elecPreSkim

if runTopSkim:
    process.pemu50_35 += process.muPreSkim
    if not runMuonSkim:
        process.plep35_35_prescale += process.boolFalse
        process.plep50_35          += process.boolFalse
        process.plep35_35_mLL170   += process.boolFalse
else:
    process.pemu50_35 += process.boolFalse


process.hnuSkimFilter = cms.EDFilter("HeavyNuFilter",
    muonTag     = cms.InputTag("patMuons"),
    electronTag = cms.InputTag("patElectrons"),
    jetTag      = cms.InputTag("selectedPatJetsPFlow"), 
    electronRho = cms.InputTag("kt6PFJetsForIsolation","rho"),

    mode = cms.int32( 991 ),
    minNmuons = cms.int32( 2 ), 
    minNelecs = cms.int32( 0 ), 

    minMu1pt = cms.double(35),
    minMu2pt = cms.double(35),
    minJetPt = cms.double(0),

    maxMuAbsEta  = cms.double(2.4),
    maxElecAbsEta = cms.double(2.5),
    maxJetAbsEta = cms.double(2.5),
    minMuMuMass  = cms.double(0),
    min4objMass  = cms.double(0),
    minMuonJetdR = cms.double(0.5),
         
    isPFJets = cms.bool(True)
)
#--- Two muons --> prescaled skim
process.hnuSkimFilterMuMu = process.hnuSkimFilter.clone()
process.hnuSkimFilterMuMu.minNmuons = cms.int32(2)
process.hnuSkimFilterMuMu.minNelecs = cms.int32(0)
process.hnuSkimFilterMuMu.mode      = cms.int32(991)
process.hnuSkimFilterMuMu.minMu1pt  = cms.double(35)
process.hnuSkimFilterMuMu.minMu2pt  = cms.double(35)
#--- Two muons, high pT
process.hnuSkimFilterMuMuPt = process.hnuSkimFilter.clone()
process.hnuSkimFilterMuMuPt.minNmuons = cms.int32(2)
process.hnuSkimFilterMuMuPt.minNelecs = cms.int32(0)
process.hnuSkimFilterMuMuPt.mode      = cms.int32(991)
process.hnuSkimFilterMuMuPt.minMu1pt  = cms.double(50)
process.hnuSkimFilterMuMuPt.minMu2pt  = cms.double(35)
#--- Two muons, two jets
process.hnuSkimFilterMuMuJJ = process.hnuSkimFilter.clone()
process.hnuSkimFilterMuMuJJ.minNmuons = cms.int32(2)
process.hnuSkimFilterMuMuJJ.minNelecs = cms.int32(0)
process.hnuSkimFilterMuMuJJ.mode      = cms.int32(992)
process.hnuSkimFilterMuMuJJ.minMu1pt  = cms.double(50)
process.hnuSkimFilterMuMuJJ.minMu2pt  = cms.double(35)
process.hnuSkimFilterMuMuJJ.minJetPt  = cms.double(25)
#--- Two muons, m(MuMu) > 170
process.hnuSkimFilterMuMuMass = process.hnuSkimFilter.clone() 
process.hnuSkimFilterMuMuMass.minNmuons   = cms.int32(2)
process.hnuSkimFilterMuMuMass.minNelecs   = cms.int32(0)
process.hnuSkimFilterMuMuMass.mode        = cms.int32(993)
process.hnuSkimFilterMuMuMass.minMuMuMass = cms.double(170)

#--- Two electrons --> prescaled skim
process.hnuSkimFilterEE = process.hnuSkimFilter.clone()
process.hnuSkimFilterEE.minNmuons = cms.int32(0)
process.hnuSkimFilterEE.minNelecs = cms.int32(2)
process.hnuSkimFilterEE.mode      = cms.int32(991)
process.hnuSkimFilterEE.minMu1pt  = cms.double(35)
process.hnuSkimFilterEE.minMu2pt  = cms.double(35)
#--- Two electrons, high pT 
process.hnuSkimFilterEEpt = process.hnuSkimFilter.clone()
process.hnuSkimFilterEEpt.minNmuons = cms.int32(0)
process.hnuSkimFilterEEpt.minNelecs = cms.int32(2)
process.hnuSkimFilterEEpt.mode      = cms.int32(991)
process.hnuSkimFilterEEpt.minMu1pt  = cms.double(50)
process.hnuSkimFilterEEpt.minMu2pt  = cms.double(35)
#--- Two electrons, two jets
process.hnuSkimFilterEEJJ = process.hnuSkimFilter.clone()
process.hnuSkimFilterEEJJ.minNmuons = cms.int32(0)
process.hnuSkimFilterEEJJ.minNelecs = cms.int32(2)
process.hnuSkimFilterEEJJ.mode      = cms.int32(992)
process.hnuSkimFilterEEJJ.minMu1pt  = cms.double(50)
process.hnuSkimFilterEEJJ.minMu2pt  = cms.double(35)
process.hnuSkimFilterEEJJ.minJetPt  = cms.double(25)
#--- Two electrons, m(EE) > 170
process.hnuSkimFilterEEMass = process.hnuSkimFilter.clone() 
process.hnuSkimFilterEEMass.minNmuons   = cms.int32(0)
process.hnuSkimFilterEEMass.minNelecs   = cms.int32(2)
process.hnuSkimFilterEEMass.mode        = cms.int32(993)
process.hnuSkimFilterEEMass.minMuMuMass = cms.double(170)

#--- Top-like: One muon, one electron, two jets
process.hnuSkimFilterTop = process.hnuSkimFilter.clone()
process.hnuSkimFilterTop.minNmuons = cms.int32(1)
process.hnuSkimFilterTop.minNelecs = cms.int32(1)
process.hnuSkimFilterTop.mode      = cms.int32(992)
process.hnuSkimFilterTop.minMu1pt  = cms.double(50)
process.hnuSkimFilterTop.minMu2pt  = cms.double(35)
process.hnuSkimFilterTop.minJetPt  = cms.double(25)

process.prescaleFilter = cms.EDFilter("Prescaler",
    prescaleFactor = cms.int32(thePrescaleFactor),
    prescaleOffset = cms.int32(0)
)

if isPFJets:
    process.plep35_35_prescale += process.modifiedPF2PATSequence
    process.plep50_35          += process.modifiedPF2PATSequence
    process.plep35_35_mLL170   += process.modifiedPF2PATSequence
    process.pemu50_35          += process.modifiedPF2PATSequence

if runMuonSkim:
    process.plep35_35_prescale += process.hnuSkimFilterMuMu
    process.plep50_35          += process.hnuSkimFilterMuMuPt
    process.plep35_35_mLL170   += process.hnuSkimFilterMuMuMass
if runElectronSkim:
    process.plep35_35_prescale += process.hnuSkimFilterEE
    process.plep50_35          += process.hnuSkimFilterEEpt
    process.plep35_35_mLL170   += process.hnuSkimFilterEEMass
    
if runTopSkim:
    process.pemu50_35 += process.hnuSkimFilterTop
    
process.plep35_35_prescale += process.prescaleFilter

if isData:
    from PhysicsTools.PatAlgos.tools.coreTools import *
    if isPFJets:
        removeMCMatchingPF2PAT( process, '' )
    if cmsswRelease == 53:
        removeMCMatching(process, ['All'], outputModules = [])    

## ============================== ##
## Python tools --> Order matters ##
## ============================== ##
from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
from PhysicsTools.PatAlgos.tools.coreTools import removeSpecificPATObjects
if cmsswRelease == 53:
    removeCleaning( process, outputModules = [] )
    removeSpecificPATObjects(process, names = ['Jets','Taus','METs'], outputModules = [])

#--- Restore the event content after PAT ---#
process.out.outputCommands = cms.untracked.vstring("keep *","drop *_*_*_HNUSKIMS","keep edmTriggerResults_TriggerResults__HNUSKIMS")

#----------------#
#--- LumiList ---#
#----------------#
process.lumilist = cms.EDAnalyzer('LumiList')
if not isMC:
    process.llPath = cms.Path(process.lumilist + process.boolFalse)

process.end = cms.EndPath(process.out)

