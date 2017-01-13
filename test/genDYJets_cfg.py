import FWCore.ParameterSet.Config as cms

process = cms.Process("GenDYJetsHTandInclusive")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('ExoAnalysis.cmsWR.genDYJetsElectronChannelModules_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

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

#################################
#Paths
process.studyDYJetsDecay = cms.Path(
		#process.printDYJetsParticleTree
		process.dyJetsBareMatchedGenGluon
		*process.dyJetsBareMatchedGenQuark
		*process.dyJetsMergeGenMatchedPartons #merge gen quark and gen gluon collections
		*process.genHTFilter
		*process.genMatchedAnalyzerOne
		)

process.schedule = cms.Schedule(process.studyDYJetsDecay)


process.TFileService = cms.Service("TFileService",
		#fileName = cms.string('hadronKinematicsGenDYJetsMadHT100to200.root')
		fileName = cms.string('hadronKinematicsGenDYJetsMadIncl.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

process.source = cms.Source( "PoolSource",
	#fileNames = cms.untracked.vstring(
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAOD2015/DYMadHT100to200MiniAOD1.root',
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAOD2015/DYMadHT100to200MiniAOD2.root'
	#	),
	fileNames = cms.untracked.vstring(
		'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAOD2015/DYMadInclusiveMiniAOD1.root',
		'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAOD2015/DYMadInclusiveMiniAOD2.root',
		),
	#inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)



