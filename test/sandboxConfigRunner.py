import FWCore.ParameterSet.Config as cms

process = cms.Process("skim")

## load the filters, producers, and sequences defined in other config file fragments
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(8000)

from ExoAnalysis.cmsWR.additionalVarParsing_cff import *


#################################
#Producers
process.addStringIdentifier = cms.EDProducer('produceStringTag',
		outputCollectionName = cms.string("datasetIdentifier"),
		stringStoredInOutputCollection = cms.string(options.datasetTag)
		)

outFileName='file:'
outFileName += options.datasetTag
outFileName +='.root'
print outFileName

#################################
#Output modules
process.outputEvts = cms.OutputModule("PoolOutputModule",
		compressionAlgorithm = cms.untracked.string('LZMA'),
		compressionLevel = cms.untracked.int32(4),
		dataset = cms.untracked.PSet(
			dataTier = cms.untracked.string('MINIAODSIM'),
			filterName = cms.untracked.string('')
			),
		dropMetaData = cms.untracked.string('ALL'),
		eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
		fastCloning = cms.untracked.bool(False),
		fileName = cms.untracked.string(outFileName),
		overrideInputFileSplitLevels = cms.untracked.bool(True),
		SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('mergeDatasets'))

		)



#################################
#Paths
process.mergeDatasets = cms.Path(
		process.addStringIdentifier
		)
process.outputEvts_step = cms.EndPath(process.outputEvts)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

process.source = cms.Source( "PoolSource",
	fileNames = cms.untracked.vstring('/store/user/shervin/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJets_100to200_SHv2/160124_155404/0000/output_1.root'),
	inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


