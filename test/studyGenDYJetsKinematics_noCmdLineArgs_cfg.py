import FWCore.ParameterSet.Config as cms

process = cms.Process("StudyGenDYJetsKinematics")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('ExoAnalysis.cmsWR.genDYJetsElectronChannelModules_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

#################################
#Analyzers

## analyze the kinematic distributions of electrons and jets at gen level with these analyzers
process.genMatchedGluonAnalyzer = cms.EDAnalyzer('singleParticleAnalyzer',
		treeName = cms.string("genGluons"),
		rightHandWsCollection = cms.InputTag("dyJetsBareMatchedGenGluon")
		)

process.genMatchedQuarkAnalyzer = cms.EDAnalyzer('singleParticleAnalyzer',
		treeName = cms.string("genQuarks"),
		rightHandWsCollection = cms.InputTag("dyJetsBareMatchedGenQuark")
		)

process.genMatchedAnalyzerOne = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("genLeptonsAndJetsNoCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsCollection = cms.InputTag("dyJetsBareMatchedGenEle"), 
		jetsCollection = cms.InputTag("dyJetsMatchGenJetsToGenPartons","matchedGenJets")
		)

process.genMatchedAnalyzerTwo = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("genLeptonsAndJetsAfterPtEtaDrCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsCollection = cms.InputTag("dyJetsMatchedGenJetLeptonDrSeparation","leptonsPassingDrCut"), 
		jetsCollection = cms.InputTag("dyJetsMatchedGenJetLeptonDrSeparation","jetsPassingDrCut")
		)

process.genMatchedAnalyzerThree = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("genLeptonsAndJetsAfterAllCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsCollection = cms.InputTag("dyJetsMatchedGenFourObjMass","leptonsPassingFourObjMassCut"), 
		jetsCollection = cms.InputTag("dyJetsMatchedGenFourObjMass","jetsPassingFourObjMassCut")
		)


#################################
#Paths
process.studyDYJetsDecay = cms.Path(
		#process.printDYJetsParticleTree
		process.dyJetsBareMatchedGenLeptonsAndPartonsSeq
		*process.dyJetsMatchGenJetsToGenPartonsSeq
		*process.genMatchedAnalyzerOne
		*process.dyJetsPtEtaRestrictedMatchedGenElesAndJetsSeq
		*process.dyJetsMatchedGenJetLeptonDrSeparationSeq
		*process.genMatchedAnalyzerTwo
		*process.dyJetsRequireHighPtMatchedGenLeptonSeq
		*process.dyJetsSelectLeptonsAfterDrSeparationSeq
		*process.dyJetsMatchedGenDileptonCandSeq
		*process.dyJetsMatchedGenFourObjMassSeq
		*process.genMatchedAnalyzerThree
		)

process.studyGenGluons = cms.Path(
		process.dyJetsBareMatchedGenGluon*process.dyJetsBareMatchedGenGluonFilter*process.genMatchedGluonAnalyzer
		)

process.studyGenQuarks = cms.Path(
		process.dyJetsBareMatchedGenQuark*process.dyJetsBareMatchedGenQuarkFilter*process.genMatchedQuarkAnalyzer
		)

process.schedule = cms.Schedule(process.studyDYJetsDecay,process.studyGenGluons,process.studyGenQuarks)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('genDYJetsDecayKinematics.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

process.source = cms.Source( "PoolSource",
	fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODBkgndFiles/DYJetsToLL_M_200to400_amcnlo_13TeV_25ns_reMiniAOD_26kevts.root'),
	#inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)



