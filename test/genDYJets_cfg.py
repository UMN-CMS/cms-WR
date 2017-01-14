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

process.genMatchedAnalyzerTwo = cms.EDAnalyzer('hadronAnalyzer',
		treeName = cms.string("genElectronsNoCuts"),
		inputCollection = cms.InputTag("dyJetsBareMatchedGenEle")
		)

#################################
#Paths
process.studyDYJetsDecay = cms.Path(
		#process.printDYJetsParticleTree
		process.dyJetsBareMatchedGenGluon
		*process.dyJetsBareMatchedGenQuark
		*process.dyJetsMergeGenMatchedPartons #merge gen quark and gen gluon collections
		#*process.genHTFilter
		#select electrons from Z then calculate Z pT
		*process.dyJetsBareMatchedGenEle
		*process.dyJetsBareMatchedGenEleFilter
		*process.genMatchedAnalyzerOne
		*process.genMatchedAnalyzerTwo
		)

process.schedule = cms.Schedule(process.studyDYJetsDecay)


process.TFileService = cms.Service("TFileService",
		#fileName = cms.string('hadronAndEleKinematicsGenDYJetsMadHT600toInf.root')
		#fileName = cms.string('hadronAndEleKinematicsGenDYJetsMadHT400to600.root')
		fileName = cms.string('hadronAndEleKinematicsGenDYJetsMadHT200to400.root')
		#fileName = cms.string('hadronAndEleKinematicsGenDYJetsMadHT100to200.root')
		#fileName = cms.string('hadronAndEleKinematicsGenDYJetsMadIncl.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

process.source = cms.Source( "PoolSource",
	#fileNames = cms.untracked.vstring(
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT600toInfMiniAOD0.root',
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT600toInfMiniAOD1.root'
	#	),
	#fileNames = cms.untracked.vstring(
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT400to600MiniAOD1.root',
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT400to600MiniAOD0.root'
	#	),
	fileNames = cms.untracked.vstring(
		'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT200to400MiniAOD1.root',
		'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT200to400MiniAOD0.root'
		),
	#fileNames = cms.untracked.vstring(
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT100to200MiniAOD1.root',
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadHT100to200MiniAOD2.root'
	#	),
	#fileNames = cms.untracked.vstring(
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadInclusiveMiniAOD1.root',
	#	'file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/DYMiniAODandMinitrees2015/DYMadInclusiveMiniAOD2.root',
	#	),
	#inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)



