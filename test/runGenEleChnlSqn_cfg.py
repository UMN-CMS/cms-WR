import FWCore.ParameterSet.Config as cms

process = cms.Process("GENEEJJ")

## load the filters, producers, and sequences defined in the config file fragment
process.load('ExoAnalysis.cmsWR.linkGenSequencesForTwoDimLimits_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('standard') 
options.maxEvents = -1
options.parseArguments()

#both cutEfficiencyAndWrAndNuMasses and eeChnlCutSeq are defined in linkGenSequences cff
process.runAnalyzer = cms.Path(process.cutEfficiencyAndWrAndNuMasses)
process.runMatchedGenAnalysis = cms.Path(process.eeChnlCutSeq)

process.schedule = cms.Schedule(process.runMatchedGenAnalysis,process.runAnalyzer)

process.TFileService = cms.Service("TFileService",
		#fileName = cms.string('analysis_genElectronChannel_MWR_800_MNu_400_using_miniAOD.root')
		fileName = cms.string(options.output)

)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)


process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(options.files),
    #inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)



