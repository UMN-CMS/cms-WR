import FWCore.ParameterSet.Config as cms

process = cms.Process("CheckWRDecay")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.genElectronChannelModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('standard') 
options.maxEvents = -1
options.parseArguments()


#################################
#Analyzers

## analyze the kinematic distributions of electrons at gen and reco level with these analyzers
process.genWRAnalyzerOne = cms.EDAnalyzer('singleParticleAnalyzer',
		treeName = cms.string("genWR"),
		rightHandWsCollection = cms.InputTag("bareMatchedWR")
		)

#################################
#Paths
process.checkWRdecay = cms.Path(
		#process.egmGsfElectronIDSequence
		#*process.HEEPIDSequence
		#process.printParticleTree
		process.bareMatchedWRSeq
		*process.genWRAnalyzerOne
		#*process.selectHEEPElesSeq
		#*process.matchedRecoElesSeq
		#*process.recoAnalyzerOne
		)
process.schedule = cms.Schedule(process.checkWRdecay)


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


