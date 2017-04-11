import FWCore.ParameterSet.Config as cms

process = cms.Process('SKIM')

############################################################ OPTIONS
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('standard') 

options.register('saveTnP',
                 0, #default value = False
                 VarParsing.VarParsing.multiplicity.singleton, #singleton or list
                 VarParsing.VarParsing.varType.int, # string, int or float
                 "tell to the output module to save the tagAndProbe triggered events"
)

options.register('doWR',
		'',
		VarParsing.VarParsing.multiplicity.singleton,
		VarParsing.VarParsing.varType.string,
		"enable filters for WR skims")

options.register('GT',
		'80X_mcRun2_asymptotic_2016_miniAODv2_v1',
		VarParsing.VarParsing.multiplicity.singleton,
		VarParsing.VarParsing.varType.string,
		"global tag name")

#default options
options.maxEvents = -1

#use this option for accessing files published with CRAB
#options.files='root://xrootd-cms.infn.it/INPTFILE'

#use this option for accessing a local file
#options.files='file:WR_M-MASSWR_ToLNu_M-MASSNU_miniAOD_13TeV-2016_CHNL_NNN.root'

#options.output = "TAGNAME_skimPartNNN.root"
#

options.parseArguments()

print options


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')


process.GlobalTag.globaltag = options.GT

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.files),
                            secondaryFileNames = cms.untracked.vstring(options.secondaryFiles)
)

process.options = cms.untracked.PSet(
#    allowUnscheduled = cms.untracked.bool(False),
    wantSummary = cms.untracked.bool(True),
)

process.MessageLogger.cerr.FwkReport.reportEvery = 5000


############################################################ OUTPUT MODULES
process.load('ExoAnalysis.cmsWR.microAOD_Output_cff')
if (options.saveTnP==1):
                SelectEventsPSet = cms.untracked.PSet(
                    SelectEvents = cms.vstring( [ 'skimPreselected', 'tagAndProbe' ] )
                )
else:
                SelectEventsPSet = cms.untracked.PSet(
                    SelectEvents = cms.vstring( [ 'skimPreselected' ] )
                )


#define a process attribute for outputting a file which will be changed in a clone() call below
process.microAOD_output = cms.OutputModule("PoolOutputModule",
		compressionAlgorithm = cms.untracked.string('LZMA'),
		compressionLevel = cms.untracked.int32(4),
		dataset = cms.untracked.PSet(
			dataTier = cms.untracked.string('MINIAODSIM'),
			filterName = cms.untracked.string('')
			),
		dropMetaData = cms.untracked.string('ALL'),
		eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
		fastCloning = cms.untracked.bool(False),
		fileName = cms.untracked.string(options.output),
		outputCommands = process.MICROAODSIMEventContent.outputCommands,
		overrideInputFileSplitLevels = cms.untracked.bool(True),
		SelectEvents = SelectEventsPSet
		)


# here the full set of sequences and hlt paths used to make the first step
process.load('ExoAnalysis.cmsWR.microAOD_cff')
process.load('ExoAnalysis.cmsWR.genFilterForWrSkims_cff')

# Path and EndPath definitions
process.tagAndProbe = cms.Path(process.tagAndProbeHltSequence)
process.skimPreselected = cms.Path(process.signalHltSequence * process.wRdiLeptonSkimSequence * process.wRdijetSkimSequence)
if (options.doWR=='WRMuMu'):
	process.skimPreselected = cms.Path(process.skipGenNuEleSeq * process.skipGenNuTauSeq * process.signalHltSequence * process.wRdiLeptonSkimSequence * process.wRdijetSkimSequence)

if (options.doWR=='WREE'):
	process.skimPreselected = cms.Path(process.skipGenNuMuSeq * process.skipGenNuTauSeq * process.signalHltSequence * process.wRdiLeptonSkimSequence * process.wRdijetSkimSequence)

process.microAODoutput_step = cms.EndPath(process.microAOD_output)

############################################################ SCHEDULE
process.schedule = cms.Schedule(process.skimPreselected, process.tagAndProbe, process.microAODoutput_step)

