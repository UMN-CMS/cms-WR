import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEEJJUnmatched")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoElectronChannelUnmatchedModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#################################
#Analyzers

## analyze the kinematic distributions of the reco jets and leptons in WR->ENu->EEJJ evts with these analyzers

process.unmatchedSignalRecoAnalyzerOne = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("unmatchedSignalRecoObjectsNoCuts"),
		leptonsCollection = cms.InputTag("bareRecoLepton"),
		jetsCollection = cms.InputTag("bareRecoJet"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1)
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



#################################
#Paths
process.unmatchedRecoSignalPath = cms.Path(
		process.bareRecoParticleSeq
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



#process.schedule = cms.Schedule(process.unmatchedRecoSignalPath)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('analysis_recoElectronChannel_no_matching.root')

)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
		##WR->ENu->EEJJ 13TeV, 40 PU, MWR=2.6 TeV, MNu=1.3 TeV miniAOD files
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_1.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_23.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_17.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_24.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_16.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_13.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_22.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_12.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_20.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_26.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_10.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_11.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_18.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_25.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_14.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_21.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_4.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_19.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_15.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_9.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_6.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_3.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_2.root',

    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)


