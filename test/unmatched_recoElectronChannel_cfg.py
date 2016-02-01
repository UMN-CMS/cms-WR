import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEEJJUnmatched")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoElectronChannelSignalUnmatchedModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

from ExoAnalysis.cmsWR.additionalVarParsing_cff import *

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GT, '')

##jet energy corrections
from ExoAnalysis.cmsWR.JEC_cff import *
JEC_correction(process, options.GT)

from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")

#################################
#Filters
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.trigFilt = hltHighLevel.clone()
process.trigFilt.HLTPaths = ['HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*']
process.trigFilt.andOr = True  #if True, then multiple HLT paths will be combined with OR logic


#################################
#Analyzers

## analyze the kinematic distributions of the reco jets and leptons in WR->ENu->EEJJ evts with these analyzers

process.unmatchedSignalRecoAnalyzerOne = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("unmatchedSignalRecoObjectsNoCuts"),
		leptonsCollection = cms.InputTag("bareRecoLepton"),
		jetsCollection = cms.InputTag("bareRecoJet"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		)

process.unmatchedSignalRecoAnalyzerTwo = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("unmatchedSignalRecoObjectsWithPtEtaCuts"),
		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		jetsCollection = cms.InputTag("ptEtaRestrictedRecoJets"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1)
		)

process.unmatchedSignalRecoAnalyzerThree = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("unmatchedSignalRecoObjectsWithPtEtaAndDileptonMassCuts"),
		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		jetsCollection = cms.InputTag("ptEtaRestrictedRecoJets"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.0)
		)

process.unmatchedSignalRecoAnalyzerFour = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("testUnmatchedSignalRecoObjectsWithPtEtaDileptonMassAndDrCuts"),
		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		jetsCollection = cms.InputTag("recoDrSeparation","jetsPassingDrSeparationCut"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.0)
		)

process.unmatchedSignalRecoAnalyzerFive = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("signalRecoObjectsWithAllCuts"),
		leptonsCollection = cms.InputTag("recoFourObjMass","leadingLeptonsPassingFourObjMassCut"),
		jetsCollection = cms.InputTag("recoFourObjMass","leadingJetsPassingFourObjMassCut"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.0)
		)

#################################
#Paths
process.unmatchedRecoSignalPath = cms.Path(
		process.patJetCorrFactorsReapplyJEC
		+process.patJetsReapplyJEC
		+process.trigFilt
		*process.egmGsfElectronIDSequence
		*process.HEEPIDSequence
		*process.bareRecoParticleSeq
		*process.bareRecoDrSeparationSeq
		#*process.unmatchedSignalRecoAnalyzerOne
		*process.ptEtaRestrictedSeq
		#*process.unmatchedSignalRecoAnalyzerTwo
		*process.recoDileptonCandidateSeq
		#*process.unmatchedSignalRecoAnalyzerThree
		#*process.unmatchedSignalRecoAnalyzerFour
		*process.recoFourObjMassSeq
		*process.unmatchedSignalRecoAnalyzerFive
		)

#process.unmatchedSignalTriggerStudyPath = cms.Path(
#		process.investigateTrigger
#		)


process.schedule = cms.Schedule(process.unmatchedRecoSignalPath)


process.TFileService = cms.Service("TFileService",
		#fileName = cms.string('/eos/uscms/store/user/skalafut/DoubleEG/analyzed_DoubleEG_skim_signal_region_eejj.root')
		fileName = cms.string(options.output)
	
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

inputFiles = cms.untracked.vstring(options.files)


process.source = cms.Source( "PoolSource",
    #fileNames = cms.untracked.vstring(
    #),
	fileNames = inputFiles,
    #inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


