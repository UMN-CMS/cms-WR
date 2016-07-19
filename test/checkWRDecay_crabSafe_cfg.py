import FWCore.ParameterSet.Config as cms

process = cms.Process("CheckWRDecay")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.genElectronChannelModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10)

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
		jetsCollection = cms.InputTag("matchGenJetsToGenQuarksNoCuts","matchedGenJetsNoCuts"),
		dontCheckSelector = cms.bool(True),
		ignoreJets = cms.bool(False)
		)

process.genMatchedParticleAnalyzerTwo = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("genLeptonsAndJetsWithPtEtaCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsOneCollection = cms.InputTag("simultaneousPtEtaCutMatchedLeadingGenEle"),#lepton with WR mother 
		leptonsTwoCollection = cms.InputTag("simultaneousPtEtaCutMatchedSubleadingGenEle"),#lepton with Nu mother
		jetsCollection = cms.InputTag("simultaneousPtEtaCutMatchedGenJets"),
		dontCheckSelector = cms.bool(True),
		ignoreJets = cms.bool(False)
		)

process.genMatchedParticleAnalyzerTwoPFive = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("genLeptonsAndJetsWithPtEtaDrCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsOneCollection = cms.InputTag("genMatchedJetLeptonDrSeparation","matchedLeadingGenElePassingDrSeparationCut"),#lepton with WR mother 
		leptonsTwoCollection = cms.InputTag("genMatchedJetLeptonDrSeparation","matchedSubleadingGenElePassingDrSeparationCut"),#lepton with Nu mother
		jetsCollection = cms.InputTag("genMatchedJetLeptonDrSeparation","genJetsPassingDrSeparationCut"),
		dontCheckSelector = cms.bool(True),
		ignoreJets = cms.bool(False)
		)

process.genMatchedParticleAnalyzerThree = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("genLeptonsAndJetsWithPtEtaDrHighPtCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsOneCollection = cms.InputTag("pickGenMatchedLeadingEle"),#lepton with WR mother 
		leptonsTwoCollection = cms.InputTag("pickGenMatchedSubleadingEle"),#lepton with Nu mother
		jetsCollection = cms.InputTag("genMatchedJetLeptonDrSeparation","genJetsPassingDrSeparationCut"),
		dontCheckSelector = cms.bool(True),
		ignoreJets = cms.bool(False)
		)

process.genMatchedParticleAnalyzerFour = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("genLeptonsAndJetsWithPtEtaDrHighPtDileptonMassCuts"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200),
		leptonsOneCollection = cms.InputTag("pickGenMatchedLeadingEle"),#lepton with WR mother 
		leptonsTwoCollection = cms.InputTag("pickGenMatchedSubleadingEle"),#lepton with Nu mother
		jetsCollection = cms.InputTag("genMatchedJetLeptonDrSeparation","genJetsPassingDrSeparationCut"),
		dontCheckSelector = cms.bool(True),
		ignoreJets = cms.bool(False)
		)

process.genMatchedParticleAnalyzerFive = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("genLeptonsAndJetsWithAllCuts"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200),
		leptonsOneCollection = cms.InputTag("genMatchedFourObjMass","genMatchedLeadingElectronPassingFourObjMassCut"),#lepton with WR mother 
		leptonsTwoCollection = cms.InputTag("genMatchedFourObjMass","genMatchedSubleadingElectronPassingFourObjMassCut"),#lepton with Nu mother
		jetsCollection = cms.InputTag("genMatchedFourObjMass","genJetsPassingFourObjMassCut"),
		dontCheckSelector = cms.bool(True),
		ignoreJets = cms.bool(False)
		)


#################################
#Paths
process.checkWRdecay = cms.Path(
		#process.hasGenMuOrTauFlavorsSeq
		process.bareMatchedWRSeq
		*process.genWRAnalyzerOne
		*process.bareMatchedNuSeq
		*process.genNuAnalyzerOne
		##identify the gen leptons from the WR decay, and the gen jets
		##matched to the gen quarks from the WR decay
		*process.bareMatchedGenParticleSeq
		*process.bareGenJetSeq
		*process.matchGenJetsToGenQuarksSeq
		##now that the gen leptons and gen jets have been selected
		##run an analyzer to study their kinematics before any cuts
		*process.genMatchedParticleAnalyzerOne
		##now apply pt, eta, and dR(lepton, jet) cuts to the gen leptons and gen jets
		*process.simultaneousPtEtaCutMatchedObjectsSeq
		*process.genMatchedParticleAnalyzerTwo
		*process.genMatchedJetLeptonDrSeparationSeq
		*process.genMatchedParticleAnalyzerTwoPFive
		##now apply the one lepton pt>60, dilepton mass, and four object mass cuts
		*process.pickGenMatchedEleSeq
		*process.requireGenMatchedHighPtEleSeq
		*process.genMatchedParticleAnalyzerThree
		*process.genMatchedDiLeptonCandidateSeq
		*process.genMatchedParticleAnalyzerFour
		*process.genMatchedFourObjMassSeq
		*process.genMatchedParticleAnalyzerFive	
		)
process.schedule = cms.Schedule(process.checkWRdecay)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('WR_decay_kinematics_MWR_5000_MNu_2500_8TeV.root')

)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

process.source = cms.Source( "PoolSource",
	fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/privateWRGen/WR_MWR_5000_ToENu_MNu_2500_GEN_8TeV_1.root'),
	
	#inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)


