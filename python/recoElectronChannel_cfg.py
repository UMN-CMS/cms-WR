import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEEJJ")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.genElectronChannelModules_cff')
process.load('ExoAnalysis.cmsWR.recoElectronChannelModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

## analyze the kinematic distributions of the reco jets and leptons before any kinematic or dR cuts are applied
process.matchedRecoAnalyzerOne = cms.EDAnalyzer('genMatchedAnalyzer',
		treeName = cms.string("matchedRecoObjectsNoCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		genLeadingLeptonCollection = cms.InputTag("matchRecoEleToLeadingGenEleNoCuts","matchedLeadingRecoEleNoCuts"),
		genSubleadingLeptonCollection = cms.InputTag("matchRecoEleToSubleadingGenEleNoCuts","matchedSubleadingRecoEleNoCuts"),
		genQuarkCollection = cms.InputTag("matchRecoJetsToGenJetsNoCuts","matchedRecoJetsNoCuts")
		)

process.matchedRecoAnalyzerTwo = cms.EDAnalyzer('genMatchedAnalyzer',
		treeName = cms.string("matchedRecoObjectsWithPtEtaCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
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
		* process.matchRecoNoCutsSeq
		* process.matchedRecoAnalyzerOne
		* process.ptEtaRestrictedMatchedRecoSeq
		* process.matchedRecoAnalyzerTwo

		)

process.TFileService = cms.Service("TFileService",
		fileName = cms.string('analysis_recoElectronChannel.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_1.root',
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/miniAOD_WR_signal/WR_signal_miniAODFile_2.root',

    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


