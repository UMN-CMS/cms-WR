# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: miniAOD --filein file:rec_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco_RECO_PU.root --fileout file:EXO-Phys14DR-00009.root --mc --eventcontent MINIAODSIM --runUnscheduled --datatier MINIAODSIM --conditions auto:run2_mc --step PAT
import FWCore.ParameterSet.Config as cms

process = cms.Process('SKIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
		fileNames = cms.untracked.vstring('file:noFile.root'),
		secondaryFileNames = cms.untracked.vstring()
		)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(False),
    wantSummary = cms.untracked.bool(True),
	SkipEvent = cms.untracked.vstring('ProductNotFound')
	)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('miniAOD nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.load('ExoAnalysis.cmsWR.microAOD_Output_cff')

#define a process attribute for outputting a file which will be changed in a clone() call below
process.MINIAODSIM_signal_output = cms.OutputModule("PoolOutputModule",
		compressionAlgorithm = cms.untracked.string('LZMA'),
		compressionLevel = cms.untracked.int32(4),
		dataset = cms.untracked.PSet(
			dataTier = cms.untracked.string('MINIAODSIM'),
			filterName = cms.untracked.string('')
			),
		dropMetaData = cms.untracked.string('ALL'),
		eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
		fastCloning = cms.untracked.bool(False),
		fileName = cms.untracked.string('file:realData_emuSignalRegionSkim.root'),
		outputCommands = process.MICROAODSIMEventContent.outputCommands,
		overrideInputFileSplitLevels = cms.untracked.bool(True),
		SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eleMuSignalSkim'))

		)

#process.MINIAODSIM_sideband_output = cms.OutputModule("PoolOutputModule",
#    compressionAlgorithm = cms.untracked.string('LZMA'),
#    compressionLevel = cms.untracked.int32(4),
#    dataset = cms.untracked.PSet(
#        dataTier = cms.untracked.string('MINIAODSIM'),
#        filterName = cms.untracked.string('')
#    ),
#    dropMetaData = cms.untracked.string('ALL'),
#    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
#    fastCloning = cms.untracked.bool(False),
#    fileName = cms.untracked.string('file:skim_sideband.root'),
#    outputCommands = process.MICROAODSIMEventContent.outputCommands,
#    overrideInputFileSplitLevels = cms.untracked.bool(True),
#    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('diMuonSidebandSkim','diElectronSidebandSkim'))
#)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.load('ExoAnalysis.cmsWR.microAOD_cff')
process.load('ExoAnalysis.cmsWR.skimEMu_cff')


# Path and EndPath definitions
process.eleMuSignalSkim = cms.Path(process.emuwRdiLeptonAndFourObjSignalSeq)

#process.diMuonSidebandSkim = cms.Path(process.wRdiMuonSidebandSeq)
#process.diElectronSidebandSkim = cms.Path(process.wRdiElectronSidebandSeq)

#process.MINIAODSIMoutput_step = cms.EndPath(process.microAODslimmingSeq * (process.MINIAODSIM_signal_output + process.MINIAODSIM_sideband_output))

process.MINIAODSIMoutput_step = cms.EndPath(process.microAODslimmingSeq * process.MINIAODSIM_signal_output )


#do not add changes to your config after this point (unless you know what you are doing)

# End of customisation functions
