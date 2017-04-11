# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: --mc --step DIGI,L1,DIGI2RAW,HLT:25ns10e33_v2 --filein=file:WR-2400_ToLNu-1200_ToEEJJ_NEW-GEN-SIM.root --fileout=file:WR_M-2400_ToLNu_M-1200_HLT_NEW.root --python_filename=WR_M-2400_ToLNu_M-1200_HLT_NEW.py --eventcontent=RAWSIM --datatier=GEN-SIM-RAW --conditions=80X_mcRun2_asymptotic_2016_miniAODv2_v1 --era=Run2_2016 --pileup_input=dbs:/MinBias_TuneCUETP8M1_13TeV-pythia8/RunIISummer15GS-MCRUN2_71_V1-v2/GEN-SIM --pileup=Flat_10_50_25ns -n 100 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('HLT',eras.Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_Flat_10_50_25ns_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_25ns10e33_v2_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
		input = cms.untracked.int32(-1)
		)

# Input source
process.source = cms.Source("PoolSource",
		dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
		fileNames = cms.untracked.vstring('file:/tmp/skalafut/WR_M-MASSWR_ToLNu_M-MASSNU_GEN-SIM_13TeV-2016_EEJJ_NUM.root'),
		inputCommands = cms.untracked.vstring('keep *', 
			'drop *_genParticles_*_*', 
			'drop *_genParticlesForJets_*_*', 
			'drop *_kt4GenJets_*_*', 
			'drop *_kt6GenJets_*_*', 
			'drop *_iterativeCone5GenJets_*_*', 
			'drop *_ak4GenJets_*_*', 
			'drop *_ak7GenJets_*_*', 
			'drop *_ak8GenJets_*_*', 
			'drop *_ak4GenJetsNoNu_*_*', 
			'drop *_ak8GenJetsNoNu_*_*', 
			'drop *_genCandidatesForMET_*_*', 
			'drop *_genParticlesForMETAllVisible_*_*', 
			'drop *_genMetCalo_*_*', 
			'drop *_genMetCaloAndNonPrompt_*_*', 
			'drop *_genMetTrue_*_*', 
			'drop *_genMetIC5GenJs_*_*'),
		secondaryFileNames = cms.untracked.vstring()
		)

process.options = cms.untracked.PSet(

		)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
		annotation = cms.untracked.string('--mc nevts:100'),
		name = cms.untracked.string('Applications'),
		version = cms.untracked.string('$Revision: 1.19 $')
		)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
		dataset = cms.untracked.PSet(
			dataTier = cms.untracked.string('GEN-SIM-RAW'),
			filterName = cms.untracked.string('')
			),
		eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
		fileName = cms.untracked.string('file:WR_M-MASSWR_ToLNu_M-MASSNU_RAWSIM_13TeV-2016_EEJJ_NUM.root'),
		outputCommands = process.RAWSIMEventContent.outputCommands,
		splitLevel = cms.untracked.int32(0)
		)

# Additional output definition

# Other statements
process.mix.input.fileNames = cms.untracked.vstring(['/store/mc/RunIISummer15GS/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM/MCRUN2_71_V1-v2/10000/004CC894-4877-E511-A11E-0025905C3DF8.root', '/store/mc/RunIISummer15GS/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM/MCRUN2_71_V1-v2/10000/0063EDE9-2F77-E511-BAF6-0002C90B7F30.root', '/store/mc/RunIISummer15GS/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM/MCRUN2_71_V1-v2/10000/0091527A-3E77-E511-B123-002590AC4BF6.root', '/store/mc/RunIISummer15GS/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM/MCRUN2_71_V1-v2/10000/00BA861E-7779-E511-85DC-0024E85A3F69.root', '/store/mc/RunIISummer15GS/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM/MCRUN2_71_V1-v2/10000/00F372BD-3C77-E511-8D36-0025901E4F3C.root', '/store/mc/RunIISummer15GS/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM/MCRUN2_71_V1-v2/10000/00FCB56F-4377-E511-8F47-0025905C2CBC.root', '/store/mc/RunIISummer15GS/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM/MCRUN2_71_V1-v2/10000/02310BE5-8F79-E511-AD22-02163E010E73.root', '/store/mc/RunIISummer15GS/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM/MCRUN2_71_V1-v2/10000/023B5EF1-4177-E511-A3E7-00266CFFC7CC.root', '/store/mc/RunIISummer15GS/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM/MCRUN2_71_V1-v2/10000/02469931-4377-E511-8A79-00259048AC98.root', '/store/mc/RunIISummer15GS/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM/MCRUN2_71_V1-v2/10000/0275943C-5477-E511-A9C5-002481D24972.root'])
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_miniAODv2_v1', '')  

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RAWSIMoutput_step])

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforFullSim 

#call to customisation function customizeHLTforFullSim imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforFullSim(process)
