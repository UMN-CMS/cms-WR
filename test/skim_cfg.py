# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: miniAOD --filein file:rec_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco_RECO_PU.root --fileout file:EXO-Phys14DR-00009.root --mc --eventcontent MINIAODSIM --runUnscheduled --datatier MINIAODSIM --conditions auto:run2_mc --step PAT
import os
import FWCore.ParameterSet.Config as cms

doMuonChannel = False 
doTTBar = False 
doDyPlusJets = True

print 'doTTBar is ', doTTBar
print 'doDyPlusJets is ', doDyPlusJets
print 'doMuonChannel is ', doMuonChannel 

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

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.load('ExoAnalysis.cmsWR.skim_cff')

# Input sources for muon and ele channels, both signal and bkgnd
if(doMuonChannel == False and doTTBar == False and doDyPlusJets == False):
	process.source = process.source.clone(
			fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuEToEEJJ_MW-800_MNu-400_TuneCUETP8M1_pythia8_13TeV_1.root')		
			)##end clone()

#

if(doMuonChannel == False and doTTBar == True):
	process.source = process.source.clone(
			fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODBkgndFiles/TTJets_TuneCUETP8M1_13TeV_pythia8_1.root')
			)##end clone()

#

if(doMuonChannel == False and doDyPlusJets == True):
	process.source = process.source.clone(
			fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODBkgndFiles/DYJetsToLL_M-50_TuneCUETP8M1_FlatPU_10_to_50_13TeV_pythia8_1.root')
			)##end clone()

#


if(doMuonChannel == True and doTTBar == False and doDyPlusJets == False):
	process.source = process.source.clone(
			fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuMuToMuMuJJ_MW-800_MNu-400_TuneCUETP8M1_pythia8_13TeV_1.root'),
			)##end clone()

#

if(doMuonChannel == True and doTTBar == True):
	process.source = process.source.clone(
			fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODBkgndFiles/TTJets_TuneCUETP8M1_13TeV_pythia8_1.root')
			)##end clone()

#

if(doMuonChannel == True and doDyPlusJets == True):
	process.source = process.source.clone(
			fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODBkgndFiles/DYJetsToLL_M-50_TuneCUETP8M1_FlatPU_10_to_50_13TeV_pythia8_1.root')
			)##end clone()

#

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
		fileName = cms.untracked.string('file:noOutputFile.root'),
		outputCommands = process.MICROAODSIMEventContent.outputCommands,
		overrideInputFileSplitLevels = cms.untracked.bool(True),
		SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('pathName'))

		)

if(doMuonChannel == True and doTTBar == False and doDyPlusJets == False):
	process.MINIAODSIM_signal_output = process.MINIAODSIM_signal_output.clone(
			fileName = cms.untracked.string('file:skim_MuonSignal.root'),
			SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('signalMuonSkim'))
			)##end clone
#

if(doMuonChannel == True and doTTBar == True):
	process.MINIAODSIM_signal_output = process.MINIAODSIM_signal_output.clone(
			fileName = cms.untracked.string('file:skim_MuonTTBarBkgnd.root'),
			SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('bkgndMuonSkim'))
			)##end clone
#

if(doMuonChannel == True and doDyPlusJets == True):
	process.MINIAODSIM_signal_output = process.MINIAODSIM_signal_output.clone(
			fileName = cms.untracked.string('file:skim_MuonDyPlusJetsBkgnd.root'),
			SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('bkgndMuonSkim'))
			)##end clone
#

if(doMuonChannel == False and doTTBar == False and doDyPlusJets == False):
	process.MINIAODSIM_signal_output = process.MINIAODSIM_signal_output.clone(
			fileName = cms.untracked.string('file:skim_ElectronSignal.root'),
			SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('signalElectronSkim'))
			)
#

if(doMuonChannel == False and doTTBar == True):
	process.MINIAODSIM_signal_output = process.MINIAODSIM_signal_output.clone(
			fileName = cms.untracked.string('file:skim_ElectronTTBarBkgnd.root'),
			SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('bkgndElectronSkim'))
			)
#

if(doMuonChannel == False and doDyPlusJets == True):
	process.MINIAODSIM_signal_output = process.MINIAODSIM_signal_output.clone(
			fileName = cms.untracked.string('file:skim_ElectronDyPlusJetsBkgnd.root'),
			SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('bkgndElectronSkim'))
			)
#

process.MINIAODSIM_sideband_output = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('file:skim_sideband.root'),
    outputCommands = process.MICROAODSIMEventContent.outputCommands,
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('diMuonSidebandSkim','diElectronSidebandSkim','EMuSidebandSkim'))
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.load('ExoAnalysis.cmsWR.microAOD_cff')
process.load('ExoAnalysis.cmsWR.skimElectron_cff')
process.load('ExoAnalysis.cmsWR.skimMuon_cff')
process.load('ExoAnalysis.cmsWR.skimEMu_cff')
process.load('ExoAnalysis.cmsWR.genElectronChannelModules_cff')
process.load('ExoAnalysis.cmsWR.genMuonChannelModules_cff')


# Path and EndPath definitions
process.eleMuSkim = cms.Path(process.emuwRdiLeptonAndFourObjSignalSeq)

if(doMuonChannel == False and ( doTTBar == True or doDyPlusJets == True) ):
	process.bkgndElectronSkim = cms.Path(process.wRdiElectronAndFourObjSignalSeq)

#
if(doMuonChannel == True and ( doTTBar == True or doDyPlusJets == True) ):
	process.bkgndMuonSkim = cms.Path(process.wRdiMuonAndFourObjSignalSeq)

#

if(doMuonChannel == False and doTTBar == False and doDyPlusJets == False):
	process.signalElectronSkim = cms.Path(
			process.bareMatchedGenParticleSeq
			*process.etaRestrictedMatchedGenParticleSeq
			*process.wRdiElectronAndFourObjSignalSeq
			)

#

if(doMuonChannel == True and doTTBar == False and doDyPlusJets == False):
	process.signalMuonSkim = cms.Path(
			process.muBareMatchedGenParticleSeq
			*process.muEtaRestrictedMatchedGenParticleSeq
			*process.wRdiMuonAndFourObjSignalSeq
			)

#

#process.diMuonSidebandSkim = cms.Path(process.wRdiMuonSidebandSeq)
#process.diElectronSidebandSkim = cms.Path(process.wRdiElectronSidebandSeq)

#process.MINIAODSIMoutput_step = cms.EndPath(process.microAODslimmingSeq * (process.MINIAODSIM_signal_output + process.MINIAODSIM_sideband_output))

process.MINIAODSIMoutput_step = cms.EndPath(process.microAODslimmingSeq * process.MINIAODSIM_signal_output )


#do not add changes to your config after this point (unless you know what you are doing)

# End of customisation functions
