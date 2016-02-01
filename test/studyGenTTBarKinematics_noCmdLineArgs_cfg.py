import FWCore.ParameterSet.Config as cms

process = cms.Process("StudyGenTTBarKinematics")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('ExoAnalysis.cmsWR.genTTBarElectronChannelModules_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

#################################
#Analyzers

## analyze the kinematic distributions of electrons and jets at gen level with these analyzers
process.genMatchedAnalyzerOne = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("genLeptonsAndJetsNoCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsCollection = cms.InputTag("ttBarBareMatchedGenEle"), 
		jetsCollection = cms.InputTag("ttBarMatchGenJetsToGenPartons","matchedGenJets")
		)

process.genMatchedAnalyzerTwo = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("genLeptonsAndJetsAfterPtEtaDrCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsCollection = cms.InputTag("ttBarMatchedGenJetLeptonDrSeparation","leptonsPassingDrCut"), 
		jetsCollection = cms.InputTag("ttBarMatchedGenJetLeptonDrSeparation","jetsPassingDrCut")
		)

process.genMatchedAnalyzerThree = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("genLeptonsAndJetsAfterAllCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsCollection = cms.InputTag("ttBarMatchedGenFourObjMass","leptonsPassingFourObjMassCut"), 
		jetsCollection = cms.InputTag("ttBarMatchedGenFourObjMass","jetsPassingFourObjMassCut")
		)



#################################
#Paths
process.studyTTBarDecay = cms.Path(
		#process.printTTBarParticleTree
		process.ttBarBareMatchedGenLeptonsAndPartonsSeq
		*process.ttBarMatchGenJetsToGenPartonsSeq
		*process.genMatchedAnalyzerOne
		*process.ttBarPtEtaRestrictedMatchedGenElesAndJetsSeq
		*process.ttBarMatchedGenJetLeptonDrSeparationSeq
		*process.genMatchedAnalyzerTwo
		*process.ttBarRequireHighPtMatchedGenLeptonSeq
		*process.ttBarSelectLeptonsAfterDrSeparationSeq
		*process.ttBarMatchedGenDileptonCandSeq
		*process.ttBarMatchedGenFourObjMassSeq
		*process.genMatchedAnalyzerThree
		)

process.schedule = cms.Schedule(process.studyTTBarDecay)

process.TFileService = cms.Service("TFileService",
		fileName = cms.string('genTTBarDecayKinematics.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

process.source = cms.Source( "PoolSource",
	fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODBkgndFiles/TTBar_powheg_pythia8_13TeV_25ns_reMiniAOD_50kevts.root'),
	#inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)



