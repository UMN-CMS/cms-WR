import FWCore.ParameterSet.Config as cms

process = cms.Process("checkWeights")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.factAndRenormWeight_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5)


#################################
#Analyzers

#################################
#Paths
process.checkWeightPath = cms.Path(
		process.checkWeights
		)
process.schedule = cms.Schedule(process.checkWeightPath)


process.TFileService = cms.Service("TFileService",
		#fileName = cms.string(options.output)
		fileName = cms.string('file:analyzedFactAndRenormScaleWeights.root')

)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)

#inputFiles = cms.untracked.vstring(options.files)
 
process.source = cms.Source( "PoolSource",
    #fileNames = inputFiles,
	fileNames = cms.untracked.vstring('file:DYJets_amctnlo_100_1_Q2b.root'),
    #inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)


