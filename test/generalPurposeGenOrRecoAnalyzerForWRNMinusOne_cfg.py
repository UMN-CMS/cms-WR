import FWCore.ParameterSet.Config as cms

process = cms.Process("AnalyzeWRDecay")

## load the filters, producers, and sequences defined in other config file fragments
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('ExoAnalysis.cmsWR.wrGenAndRecoForNMinusOne_cff')
process.load('ExoAnalysis.cmsWR.ttBarGenAndRecoForNMinusOne_cff')
process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4'


# import the HEEP selection modules and sequences
from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('standard') 
options.maxEvents = -1

options.register('channel',
		'',
		VarParsing.VarParsing.multiplicity.singleton,
		VarParsing.VarParsing.varType.string,
		"enable EE or MuMu channel modules")

options.parseArguments()

#these collection names must be checked for correctness before execution
genLeptonCollName="tempLeptColl"
matchedGenLeptonCollName="tempMatchedLeptColl"
matchedGenQuarkCollName="tempMatchedGenQuarkColl"
matchedGenJetCollName="tempMatchedGenJetColl"
recoLeptCollName="tempRecoLeptColl"
if (options.channel=='EE'):
	genLeptonCollName = "bareGenEle"
	matchedGenLeptonCollName = "matchedGenEleFromWR"
	matchedGenQuarkCollName = "matchedGenPartonFromWREle"
	matchedGenJetCollName = "condenseMatchedGenJetWREleCollName"
	recoLeptCollName = "heepRecoEle"
#
if (options.channel=='MuMu'):
	genLeptonCollName = "bareGenMu"
	matchedGenLeptonCollName = "matchedGenMuonFromWR"
	matchedGenQuarkCollName = "matchedGenPartonFromWRMu"
	matchedGenJetCollName = "condenseMatchedGenJetWRMuCollName"
	recoLeptCollName = "highPtIdRecoMu"
#################################
#Analyzers

## analyze the kinematic distributions of leading leptons and leading hadrons at gen or reco level with these analyzers

#gen quarks and leptons without GEN mother requirements
process.analyzerOne = cms.EDAnalyzer('generalPurposeGenAndRecoAnalyzer',
		treeName = cms.string("genKinematicsUsingGenQuarksWithoutGenMotherRequirements"),
		inputParticleCollsAreGEN = cms.bool(True),
		inputHadronsAreQuarks = cms.bool(True),
		hadronCollection = cms.InputTag("bareGenQuark"),
		leptonCollection = cms.InputTag(genLeptonCollName)
		)

#gen quarks and leptons with GEN mother requirements
process.analyzerTwo = process.analyzerOne.clone(
		treeName = cms.string("genKinematicsUsingGenQuarksWithGenMotherRequirements"),
		hadronCollection = cms.InputTag(matchedGenQuarkCollName),
		leptonCollection = cms.InputTag(matchedGenLeptonCollName)
		)

#gen jets and leptons without GEN mother requirements
process.analyzerThree = process.analyzerOne.clone(
		treeName = cms.string("genKinematicsUsingGenJetsWithoutGenMotherRequirements"),
		inputHadronsAreQuarks = cms.bool(False),
		hadronCollection = cms.InputTag("bareGenJet")
		)

#gen jets and leptons with GEN mother requirements
process.analyzerFour = process.analyzerTwo.clone(
		treeName = cms.string("genKinematicsUsingGenJetsWithGenMotherRequirements"),
		inputHadronsAreQuarks = cms.bool(False),
		hadronCollection = cms.InputTag(matchedGenJetCollName),
		)

#reco jets and leptons, no matching requirements
process.analyzerFive = cms.EDAnalyzer('generalPurposeGenAndRecoAnalyzer',
		treeName = cms.string("recoKinematics"),
		inputParticleCollsAreGEN = cms.bool(False),
		inputHadronsAreQuarks = cms.bool(False),
		hadronCollection = cms.InputTag("bareRecoJetPassingId"),
		leptonCollection = cms.InputTag(recoLeptCollName)
		)

#gen quarks and leptons in events where the correct flavor GEN leptons must be present with correct mothers
process.analyzerSix = process.analyzerOne.clone(
		treeName = cms.string("genKinematicsUsingGenQuarksWithoutGenMotherRequirementsWithGenFlavorReqs")
		)

#gen jets and leptons in events where the correct flavor GEN leptons and b quarks must be present with correct mothers
process.analyzerSeven = process.analyzerThree.clone(
		treeName = cms.string("genKinematicsUsingGenJetsWithoutGenMotherRequirementsWithGenFlavorReqs")
		)

#################################
#Paths
process.analyzeWRdecayGenQuarksWithoutMatching = cms.Path(process.analyzerOne)
process.analyzeWRdecayGenQuarksWithMatching = cms.Path(process.analyzerTwo)
process.analyzeWRdecayGenJetsWithoutMatching = cms.Path(process.analyzerThree)
process.analyzeWRdecayGenJetsWithMatching = cms.Path(process.analyzerFour)
process.analyzeWRdecayReco = cms.Path(process.analyzerFive)
process.analyzeWRdecayGenQuarksWithFlavorReq = cms.Path(process.analyzerSix)
process.analyzeWRdecayGenJetsWithFlavorReq = cms.Path(process.analyzerSeven)

#all paths need to be modified based on the lepton channel
if (options.channel=='EE'):
	process.analyzeWRdecayGenQuarksWithFlavorReq = cms.Path(process.matchedGenEleAndQuarkFromWRSeq*process.genEleChnlGenQuarksNoMatchingSequence*process.analyzerSix)
	process.analyzeWRdecayGenJetsWithFlavorReq = cms.Path(process.matchedGenEleAndQuarkFromWRSeq*process.genEleChnlGenJetsNoMatchingSequence*process.analyzerSeven)
	process.analyzeWRdecayReco = cms.Path(
			process.egmGsfElectronIDSequence*process.HEEPIDSequence
			*process.recoEleJetPassingIdSeq*process.analyzerFive
			)
	process.analyzeWRdecayGenQuarksWithoutMatching = cms.Path(process.genEleChnlGenQuarksNoMatchingSequence*process.analyzerOne)
	process.analyzeWRdecayGenQuarksWithMatching = cms.Path(process.matchedGenEleAndQuarkFromWRSeq*process.analyzerTwo)
	process.analyzeWRdecayGenJetsWithoutMatching = cms.Path(process.genEleChnlGenJetsNoMatchingSequence*process.analyzerThree)
	process.analyzeWRdecayGenJetsWithMatching = cms.Path(
			process.matchedGenEleAndQuarkFromWRSeq
			*process.genJetSeq
			*process.matchedGenJetWREleSeq
			*process.condenseGenJetWREleSeq
			*process.analyzerFour)

#
if (options.channel=='MuMu'):
	process.analyzeWRdecayGenQuarksWithFlavorReq = cms.Path(process.matchedGenMuAndQuarkFromWRSeq*process.genMuChnlGenQuarksNoMatchingSequence*process.analyzerSix)
	process.analyzeWRdecayGenJetsWithFlavorReq = cms.Path(process.matchedGenMuAndQuarkFromWRSeq*process.genMuChnlGenJetsNoMatchingSequence*process.analyzerSeven)
	process.analyzeWRdecayReco = cms.Path(process.recoMuJetPassingIdSeq*process.analyzerFive)
	process.analyzeWRdecayGenQuarksWithoutMatching = cms.Path(process.genMuChnlGenQuarksNoMatchingSequence*process.analyzerOne)
	process.analyzeWRdecayGenQuarksWithMatching = cms.Path(process.matchedGenMuAndQuarkFromWRSeq*process.analyzerTwo)
	process.analyzeWRdecayGenJetsWithoutMatching = cms.Path(process.genMuChnlGenJetsNoMatchingSequence*process.analyzerThree)
	process.analyzeWRdecayGenJetsWithMatching = cms.Path(
			process.matchedGenMuAndQuarkFromWRSeq
			*process.genJetSeq
			*process.matchedGenJetWRMuSeq
			*process.condenseGenJetWRMuSeq
			*process.analyzerFour)

#

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
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


