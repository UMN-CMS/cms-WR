import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEEJJUnmatched")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoElectronChannelUnmatchedModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")

#################################
#Filters
process.trigSelector = cms.EDFilter("triggerFilter",
		checkThisHltPath = cms.string("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v2")
		)


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
		treeName = cms.string("testUnmatchedSignalRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts"),
		leptonsCollection = cms.InputTag("recoFourObjMass","leadingLeptonsPassingFourObjMassCut"),
		jetsCollection = cms.InputTag("recoFourObjMass","leadingJetsPassingFourObjMassCut"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.0)
		)

#process.investigateTrigger = cms.EDAnalyzer('triggerAnalyzer',
#		treeName = cms.string("electronChannelSignalHLT"),
#		commaSeparatedHltPaths = cms.string(),
#		trigResultsColl = cms.InputTag("TriggerResults","","HLT"),
#		trigObjectStandAloneColl = cms.InputTag("selectedPatTrigger")
#		)


#################################
#Paths
process.unmatchedRecoSignalPath = cms.Path(
		process.trigSelector
		*process.egmGsfElectronIDSequence
		*process.HEEPIDSequence
		*process.bareRecoParticleSeq
		*process.unmatchedSignalRecoAnalyzerOne
		*process.ptEtaRestrictedSeq
		*process.unmatchedSignalRecoAnalyzerTwo
		*process.recoDileptonCandidateSeq
		*process.unmatchedSignalRecoAnalyzerThree
		*process.recoDrSeparationSeq
		*process.unmatchedSignalRecoAnalyzerFour
		*process.recoFourObjMassSeq
		*process.unmatchedSignalRecoAnalyzerFive
		)

#process.unmatchedSignalTriggerStudyPath = cms.Path(
#		process.investigateTrigger
#		)


process.schedule = cms.Schedule(process.unmatchedRecoSignalPath)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('/eos/uscms/store/user/skalafut/DoubleEG/analyzed_DoubleEG_skim_signal_region_eejj.root')
	
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)


inputFiles = cms.untracked.vstring()
inputFiles.extend([
		##skims of real data! 13 TeV, 50ns  doubleEG July 2015
	    'file:/eos/uscms/store/user/skalafut/DoubleEG/realData_DoubleEG_13TeV_50ns_eejj_signalAndLowMassRegionSkim_atFNALLPC/150727_142936/0000/realData_electronSignalRegionSkim_1.root'


		###WR->ENu->EEJJ 13TeV, 40 PU, MWR=2.6 TeV, MNu=1.3 TeV miniAOD files
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_1.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_23.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_17.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_24.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_16.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_13.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_22.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_12.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_20.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_26.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_10.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_11.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_18.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_25.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_14.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_21.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_4.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_19.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_15.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_9.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_6.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_3.root',
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_2.root',

  
	])
	#end inputFiles



process.source = cms.Source( "PoolSource",
    #fileNames = cms.untracked.vstring(
    #),
	fileNames = inputFiles,
    #inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


