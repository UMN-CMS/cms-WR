import FWCore.ParameterSet.Config as cms

process = cms.Process("AnalyzeWRDecay")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.reformatGenAndRecoCollections_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('standard') 
options.maxEvents = -1
options.parseArguments()


#################################
#Analyzers

## analyze the kinematic distributions of leading leptons and leading hadrons at gen or reco level with these analyzers
process.analyzerOne = cms.EDAnalyzer('generalPurposeGenAndRecoAnalyzer',
		treeName = cms.string("genKinematicsUsingGenQuarksWithoutGenMotherRequirements"),
		inputParticleCollsAreGEN = cms.bool(True),
		inputHadronsAreQuarks = cms.bool(True),
		hadronCollection = cms.InputTag("hadronCollName"),
		leptonCollection = cms.InputTag("leptonCollName")
		)

process.analyzerTwo = process.analyzerOne.clone(
		treeName = cms.string("genKinematicsUsingGenQuarksWithGenMotherRequirements"),
		hadronCollection = cms.InputTag("hadronCollName"),
		leptonCollection = cms.InputTag("leptonCollName")
		)

process.analyzerThree = process.analyzerOne.clone(
		treeName = cms.string("genKinematicsUsingGenJetsWithoutGenMotherRequirements"),
		inputHadronsAreQuarks = cms.bool(False),
		hadronCollection = cms.InputTag("hadronCollName")
		)

process.analyzerFour = process.analyzerThree.clone(
		treeName = cms.string("genKinematicsUsingGenJetsWithGenMotherRequirements"),
		hadronCollection = cms.InputTag("hadronCollName"),
		leptonCollection = cms.InputTag("leptonCollName")
		)

process.analyzerFive = cms.EDAnalyzer('generalPurposeGenAndRecoAnalyzer',
		treeName = cms.string("recoKinematics"),
		inputParticleCollsAreGEN = cms.bool(False),
		inputHadronsAreQuarks = cms.bool(False),
		hadronCollection = cms.InputTag("hadronCollName"),
		leptonCollection = cms.InputTag("leptonCollName")
		)

#################################
#Paths
process.analyzeWRdecay = cms.Path(
		process.genAndRecoFormatterSeq
		*process.analyzerOne
		)

process.schedule = cms.Schedule(process.analyzeWRdecay)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string(options.output)
		#fileName = cms.string('analyzedWr.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

inputFiles = cms.untracked.vstring(options.files)

process.source = cms.Source( "PoolSource",
	fileNames = inputFiles,
	#fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/miniAOD/WR_signal_MC/WRToNuEToEEJJ_MW-6000_MNu-3000_TuneCUETP8M1_pythia8_13TeV.root'),
	#inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(options.maxEvents)
    input = cms.untracked.int32(-1)
)


