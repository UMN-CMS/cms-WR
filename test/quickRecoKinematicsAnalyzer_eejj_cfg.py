import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEEJJ")

#this is a copy of bkgndElectronChannel_cfg.py

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoElectronChannelUnmatchedModules_cff')
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
#as long as one of these two triggers is fired, the event is selected
process.trigSelector = cms.EDFilter("triggerFilter",
		checkThisHltPath = cms.string("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v1"),
		alsoCheckThisHltPath = cms.string("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v1"),

		)


#################################
#Analyzers

## analyze the kinematic distributions of lepton and jet candidates in bkgnd evts with these analyzers
process.recoAnalyzerOne = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("recoObjectsNoCuts"),
		leptonsCollection = cms.InputTag("bareRecoLepton"),
		jetsCollection = cms.InputTag("bareRecoJetLeptonDrSeparation","bareJetsPassingDrSeparationCut"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		)

process.recoAnalyzerTwo = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("recoObjectsWithPtEtaCuts"),
		leptonsCollection = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		jetsCollection = cms.InputTag("ptEtaRestrictedRecoJets"),
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
		process.trigSelector
		*process.egmGsfElectronIDSequence
		*process.HEEPIDSequence
		*process.bareRecoParticleSeq
		*process.bareRecoDrSeparationSeq
		*process.recoAnalyzerOne
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


