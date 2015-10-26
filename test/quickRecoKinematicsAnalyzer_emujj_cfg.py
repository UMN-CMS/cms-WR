import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEMuJJ")

#this is a copy of bkgndElectronChannel_cfg.py

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoEMuChannelSidebandUnmatchedModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

from ExoAnalysis.cmsWR.additionalVarParsing_cff import *

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GT, '')

# import the HEEP selection modules and sequences
from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")


#################################
#Filters
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.trigFilt = hltHighLevel.clone()
#process.trigFilt.HLTPaths = ['HLT_Mu45_eta2p1_v*','HLT_Mu50_v*','HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*']
process.trigFilt.HLTPaths = ['HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v*']
process.trigFilt.andOr = True  #if True, then multiple HLT paths will be combined with OR logic

#no need for trigSelector, or triggerFilter.cc
#process.trigSelector = cms.EDFilter("triggerFilter",
#		checkThisHltPath = cms.string("HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v1"),
#		alsoCheckThisHltPath = cms.string("HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v1")
#		)

#################################
#Analyzers

## analyze the kinematic distributions of lepton and jet candidates using these modules
process.recoAnalyzerOne = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("recoObjectsNoCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsOneCollection = cms.InputTag("emuBareRecoLeptonOne"),#electrons
		leptonsTwoCollection = cms.InputTag("emuBareRecoLeptonTwo"),#muons
		jetsCollection = cms.InputTag("emuBareRecoJet"),
		)

#the unmatchedAnalyzerForMixed... will always pick the two leading jets from the input collection
#this is also done in the recoEMuChannelSidebandUnmatchedModules_cff so that only the
#two leading jets are used in the four object mass cut

#use emuBareRecoJet here so that the jet multiplicity indicates the number of jets
#with pt>40, |eta| < 2.5, and pass the loose jet ID
process.recoAnalyzerTwo = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("recoObjectsWithPtEtaCuts"),
		leptonsOneCollection = cms.InputTag("emuBareRecoJetLeptonDrSeparation","leptonOnesPassingDrSeparationCut"),#electrons
		leptonsTwoCollection = cms.InputTag("emuBareRecoJetLeptonDrSeparation","leptonTwosPassingDrSeparationCut"),#muons
		jetsCollection = cms.InputTag("emuBareRecoJet"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1)
		)

#process.recoAnalyzerThree = cms.EDAnalyzer('unmatchedAnalyzer',
#		treeName = cms.string("recoObjectsWithPtEtaAndDileptonMassCuts"),
#		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
#		jetsCollection = cms.InputTag("ptEtaRestrictedRecoJets"),
#		doDileptonMassCut = cms.bool(True),
#		minDileptonMass = cms.double(200.0)
#		)
#
#process.recoAnalyzerFour = cms.EDAnalyzer('unmatchedAnalyzer',
#		treeName = cms.string("recoObjectsWithPtEtaDileptonMassAndDrCuts"),
#		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
#		jetsCollection = cms.InputTag("recoDrSeparation","jetsPassingDrSeparationCut"),
#		doDileptonMassCut = cms.bool(True),
#		minDileptonMass = cms.double(200.0)
#		)
#
#process.recoAnalyzerFive = cms.EDAnalyzer('unmatchedAnalyzer',
#		treeName = cms.string("recoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts"),
#		leptonsCollection = cms.InputTag("recoFourObjMass","leadingLeptonsPassingFourObjMassCut"),
#		jetsCollection = cms.InputTag("recoFourObjMass","leadingJetsPassingFourObjMassCut"),
#		doDileptonMassCut = cms.bool(True),
#		minDileptonMass = cms.double(200.0)
#		)


#################################
#Paths
process.unmatchedBkgndRecoPath = cms.Path(
		process.trigFilt
		*process.egmGsfElectronIDSequence
		*process.HEEPIDSidebandSequence  #only look for 1 HEEP electron
		*process.wrTunePMuProdSeq
		*process.isHighPtMuSeq
		*process.emuBareRecoParticleSeq
		*process.recoAnalyzerOne
		*process.emuBareRecoDrSeparationSeq
		*process.emuPrepLeptonsAfterDrCutSeq
		*process.emuLeadLeptonAndLowMassSeq
		*process.recoAnalyzerTwo
		#*process.recoDileptonCandidateSeq
		#*process.recoAnalyzerThree
		#*process.recoDrSeparationSeq
		#*process.recoAnalyzerFour
		#*process.recoFourObjMassSeq
		#*process.recoAnalyzerFive
		)
process.schedule = cms.Schedule(process.unmatchedBkgndRecoPath)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string(options.output)

)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

inputFiles = cms.untracked.vstring(options.files)
 
process.source = cms.Source( "PoolSource",
    fileNames = inputFiles, 
    #inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


