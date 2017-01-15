import FWCore.ParameterSet.Config as cms

process = cms.Process("GenAndRecoDYJetsHTandInclusive")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('ExoAnalysis.cmsWR.genDYJetsElectronChannelModules_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('ExoAnalysis.cmsWR.JEC_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '74X_mcRun2_asymptotic_v4', '')

#################################
#Filters
process.genHTFilter = cms.EDFilter('htFilter',
		cutThreshold = cms.double(100.),
		inputCollection = cms.InputTag("dyJetsMergeGenMatchedPartons")
		)

#################################
#Analyzers

## analyze the kinematic distributions of electrons and jets at gen level with these analyzers
process.genMatchedAnalyzerOne = cms.EDAnalyzer('hadronAnalyzer',
		treeName = cms.string("genGluonsAndQuarksNoCuts"),
		inputCollection = cms.InputTag("dyJetsMergeGenMatchedPartons")
		)

process.genMatchedAnalyzerTwo = cms.EDAnalyzer('hadronAnalyzer',
		treeName = cms.string("genElectronsNoCuts"),
		inputCollection = cms.InputTag("dyJetsBareMatchedGenEle")
		)

process.recoAnalyzerOne = cms.EDAnalyzer('hadronAnalyzer',
		treeName = cms.string("recoJetsNoCuts"),
		inputCollection = cms.InputTag("bareRecoJet")
		)


#################################
#Paths
process.studyDYJetsDecay = cms.Path(
		#process.printDYJetsParticleTree
		process.dyJetsBareMatchedGenGluon
		*process.dyJetsBareMatchedGenQuark
		*process.dyJetsMergeGenMatchedPartons #merge gen quark and gen gluon collections
		*process.genHTFilter
		##select electrons from Z then calculate Z pT
		#*process.dyJetsBareMatchedGenEle
		#*process.dyJetsBareMatchedGenEleFilter
		#*process.genMatchedAnalyzerOne
		#*process.genMatchedAnalyzerTwo
		#RECO sumPt of two leading jets
		*process.jecSequence
		*process.bareRecoJet
		*process.bareRecoJetFilter
		*process.recoAnalyzerOne
		)

process.schedule = cms.Schedule(process.studyDYJetsDecay)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('recoHadronKinematicsGenDYJetsMadInclPartNNN.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

process.source = cms.Source( "PoolSource",
	fileNames = cms.untracked.vstring('INPTFILE'),
	#inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)



