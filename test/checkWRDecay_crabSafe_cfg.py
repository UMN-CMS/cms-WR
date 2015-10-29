import FWCore.ParameterSet.Config as cms

process = cms.Process("CheckWRDecay")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.genElectronChannelModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

#import FWCore.ParameterSet.VarParsing as VarParsing
#options = VarParsing.VarParsing('standard') 
#options.maxEvents = -1
#options.parseArguments()


#################################
#Analyzers

## analyze the kinematic distributions of electrons at gen and reco level with these analyzers
process.genWRAnalyzerOne = cms.EDAnalyzer('singleParticleAnalyzer',
		treeName = cms.string("genWR"),
		rightHandWsCollection = cms.InputTag("bareMatchedWR")
		)

process.genNuAnalyzerOne = cms.EDAnalyzer('singleParticleAnalyzer',
		treeName = cms.string("genNu"),
		rightHandWsCollection = cms.InputTag("bareMatchedNu")
		)

process.genMatchedParticleAnalyzerOne = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("genLeptonsAndJetsNoCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsOneCollection = cms.InputTag("bareMatchedLeadingGenEle"),#lepton with WR mother 
		leptonsTwoCollection = cms.InputTag("bareMatchedSubleadingGenEle"),#lepton with Nu mother
		jetsCollection = cms.InputTag("matchGenJetsToGenQuarksNoCuts","matchedGenJetsNoCuts")
		)

process.genMatchedParticleAnalyzerTwo = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("genLeptonsAndJetsWithPtEtaCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsOneCollection = cms.InputTag("simultaneousPtEtaCutMatchedLeadingGenEle"),#lepton with WR mother 
		leptonsTwoCollection = cms.InputTag("simultaneousPtEtaCutMatchedSubleadingGenEle"),#lepton with Nu mother
		jetsCollection = cms.InputTag("simultaneousPtEtaCutMatchedGenJets")
		)

process.genMatchedParticleAnalyzerThree = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("genLeptonsAndJetsWithAllCuts"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200),
		leptonsOneCollection = cms.InputTag("genMatchedFourObjMass","genMatchedLeadingElectronPassingFourObjMassCut"),#lepton with WR mother 
		leptonsTwoCollection = cms.InputTag("genMatchedFourObjMass","genMatchedSubleadingElectronPassingFourObjMassCut"),#lepton with Nu mother
		jetsCollection = cms.InputTag("genMatchedFourObjMass","genJetsPassingFourObjMassCut")
		)


#################################
#Paths
process.checkWRdecay = cms.Path(
		process.bareMatchedWRSeq
		*process.genWRAnalyzerOne
		*process.bareMatchedNuSeq
		*process.genNuAnalyzerOne
		#identify the gen leptons from the WR decay, and the gen jets
		#matched to the gen quarks from the WR decay
		*process.bareMatchedGenParticleSeq
		*process.bareGenJetSeq
		*process.matchGenJetsToGenQuarksSeq
		#now that the gen leptons and gen jets have been selected
		#run an analyzer to study their kinematics before any cuts
		*process.genMatchedParticleAnalyzerOne
		#now apply pt and eta cuts to the gen leptons and gen jets
		*process.simultaneousPtEtaCutMatchedObjectsSeq
		*process.genMatchedParticleAnalyzerTwo
		#now apply the dR(lepton, jet), one lepton pt>60, dilepton mass, and four object mass cuts
		*process.genMatchedJetLeptonDrSeparationSeq
		*process.pickGenMatchedEleSeq
		*process.requireGenMatchedHighPtEleSeq
		*process.genMatchedDiLeptonCandidateSeq
		*process.genMatchedFourObjMassSeq
		*process.genMatchedParticleAnalyzerThree	
		)
process.schedule = cms.Schedule(process.checkWRdecay)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('genWrKinematics.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)


process.source = cms.Source( "PoolSource",
	fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuEToEEJJ_MW-1400_MNu-700-TuneCUETP8M1_pythia8_13TeV_1.root'),
	#inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


