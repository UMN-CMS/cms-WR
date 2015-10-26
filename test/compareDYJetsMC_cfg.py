import FWCore.ParameterSet.Config as cms

process = cms.Process("DYJetsComparo")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.inspectGenAndRecoZeeInDYJetsSamples_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

from ExoAnalysis.cmsWR.additionalVarParsing_cff import *

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GT, '')

from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")

#################################
#Analyzers

## analyze the kinematic distributions of electrons at gen and reco level with these analyzers
process.genAnalyzerOne = cms.EDAnalyzer('zeeAnalyzer',
		treeName = cms.string("genZedEleEle"),
		leptonsCollection = cms.InputTag("genEleFromZ"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		)

process.recoAnalyzerOne = cms.EDAnalyzer('zeeAnalyzer',
		treeName = cms.string("recoZedEleEle"),
		leptonsCollection = cms.InputTag("matchedRecoElesToGenZeeEles","matchedRecoElesFromZee"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		)


#################################
#Paths
process.compareDYJetsPath = cms.Path(
		process.egmGsfElectronIDSequence
		*process.HEEPIDSequence
		*process.genZeeSeq
		*process.genAnalyzerOne
		*process.selectHEEPElesSeq
		*process.matchedRecoElesSeq
		*process.recoAnalyzerOne
		)
process.schedule = cms.Schedule(process.compareDYJetsPath)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string(options.output)
		#fileName = cms.string('file:analyzedZToEleEleSkim.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

inputFiles = cms.untracked.vstring(options.files)

process.source = cms.Source( "PoolSource",
	fileNames = inputFiles,
	#fileNames = cms.untracked.vstring('file:noFile.root'),
	#inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    #input = cms.untracked.int32(-1)
)


