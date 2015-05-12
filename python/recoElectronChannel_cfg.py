import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEEJJ")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.genElectronChannelModules_cff')
process.load('ExoAnalysis.cmsWR.recoElectronChannelModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

## analyze the kinematic distributions of the reco jets and leptons before any kinematic or dR cuts are applied
process.matchedRecoAnalyzerOne = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedRecoObjectsNoCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		genLeadingLeptonCollection = cms.InputTag("matchRecoEleToLeadingGenEleNoCuts","matchedLeadingRecoEleNoCuts"),
		genSubleadingLeptonCollection = cms.InputTag("matchRecoEleToSubleadingGenEleNoCuts","matchedSubleadingRecoEleNoCuts"),
		genQuarkCollection = cms.InputTag("matchRecoJetsToGenJetsNoCuts","matchedRecoJetsNoCuts")
		)

process.matchedRecoAnalyzerTwo = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedRecoObjectsWithPtEtaCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		genLeadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingRecoEle"),
		genSubleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingRecoEle"),
		genQuarkCollection = cms.InputTag("ptEtaRestrictedMatchedRecoJets")
		)

process.matchedRecoAnalyzerThree = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedRecoObjectsWithPtEtaAndDileptonMassCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		genLeadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingRecoEle"),
		genSubleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingRecoEle"),
		genQuarkCollection = cms.InputTag("ptEtaRestrictedMatchedRecoJets")
		)


process.matchedRecoAnalyzerFour = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedRecoObjectsWithPtEtaDileptonMassAndDrCuts"),
		doDeltaRcut = cms.bool(True),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(0.4),
		genLeadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingRecoEle"),
		genSubleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingRecoEle"),
		genQuarkCollection = cms.InputTag("ptEtaRestrictedMatchedRecoJets")
		)


process.matchedRecoAnalyzerFive = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts"),
		doDeltaRcut = cms.bool(True),
		doFourObjMassCut = cms.bool(True),
		minFourObjMass = cms.double(600),
		minDeltaRforLeptonJetExclusion = cms.double(0.4),
		genLeadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingRecoEle"),
		genSubleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingRecoEle"),
		genQuarkCollection = cms.InputTag("ptEtaRestrictedMatchedRecoJets")
		)



process.matchedRecoPath = cms.Path(
		process.bareGenJet
		* process.bareGenJetFilter
		* process.bareMatchedGenParticleSeq
		* process.matchGenJetsToGenQuarksNoCuts
		* process.matchGenJetsToGenQuarksNoCutsFilter
		* process.bareRecoParticleSeq
		* process.testCutProducer
		* process.testCutProducerFilter
		* process.matchRecoNoCutsSeq
		* process.matchedRecoAnalyzerOne
		* process.ptEtaRestrictedMatchedRecoSeq
		* process.matchedRecoAnalyzerTwo
		* process.recoMatchedDiElectronCandidateSeq
		* process.matchedRecoAnalyzerThree
		* process.matchedRecoAnalyzerFour
		* process.matchedRecoAnalyzerFive

		)

process.onlyMatchJetsPath = cms.Path(
		process.bareGenJet
		* process.bareGenJetFilter
		* process.bareMatchedGenParticleSeq
		* process.matchGenJetsToGenQuarksNoCutsNewPath
		* process.matchGenJetsToGenQuarksNoCutsNewPathFilter
		* process.bareRecoJet
		* process.bareRecoJetFilter
		* process.matchRecoJetsToGenJetsNoCutsNewPath
		* process.matchRecoJetsToGenJetsNoCutsNewPathFilter
	
		)

process.onlyMatchLeadingElePath = cms.Path(
		process.bareMatchedGenParticleSeq
		* process.bareRecoEle
		* process.bareRecoEleFilter
		* process.matchRecoEleToLeadingGenEleNoCutsNewPath
		* process.matchRecoEleToLeadingGenEleNoCutsNewPathFilter

		)

process.onlyMatchSubleadingElePath = cms.Path(
		process.bareMatchedGenParticleSeq
		* process.bareRecoEle
		* process.bareRecoEleFilter
		* process.matchRecoEleToSubleadingGenEleNoCutsNewPath
		* process.matchRecoEleToSubleadingGenEleNoCutsNewPathFilter
		
		)



process.TFileService = cms.Service("TFileService",
		fileName = cms.string('analysis_recoElectronChannel_two_stage_matching_for_jets.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
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
    input = cms.untracked.int32(-1)
)


