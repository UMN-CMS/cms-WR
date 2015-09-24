import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEEJJ")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoElectronChannelSidebandUnmatchedModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('standard') 
options.maxEvents = -1
options.parseArguments()

from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")


#################################
#Filters
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.trigFilt = hltHighLevel.clone()
process.trigFilt.HLTPaths = ['HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*']
process.trigFilt.andOr = True  #if True, then multiple HLT paths will be combined with OR logic

#no need for trigSelector, or triggerFilter.cc
#process.trigSelector = cms.EDFilter("triggerFilter",
#		checkThisHltPath = cms.string("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v1"),
#		alsoCheckThisHltPath = cms.string("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v1"),
#
#		)


#################################
#Analyzers

## analyze the kinematic distributions of lepton and jet candidates in bkgnd evts with these analyzers
process.recoAnalyzerOne = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("recoObjectsNoCuts"),
		jetsCollection = cms.InputTag("bareRecoJet"),
		leptonsCollection = cms.InputTag("bareRecoLepton"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		)

#the unmatchedAnalyzer.cc will always pick the two leading jets from the input collection
#this is also done in the recoElectronChannelSidebandUnmatchedModules_cff so that only the
#two leading jets are used in the four object mass cut
#use bareRecoJet here so that the jet multiplicity indicates the number of jets
#with pt>40, |eta| < 2.5, and pass the loose jet ID
process.recoAnalyzerTwo = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("recoObjectsWithPtEtaCuts"),
		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		jetsCollection = cms.InputTag("bareRecoJet"),
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
		*process.HEEPIDSequence
		*process.bareRecoParticleSeq
		*process.recoAnalyzerOne
		*process.bareRecoDrSeparationSeq
		*process.ptEtaRestrictedSeq
		*process.lowMassLLJJObjectSeq
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
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


