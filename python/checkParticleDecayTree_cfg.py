import FWCore.ParameterSet.Config as cms

process = cms.Process("GENLLQQ")

## load the filters, producers, and sequences defined in the config file fragment
process.load('ExoAnalysis.cmsWR.genElectronChannelModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


process.printMuOrTauDecayChain = cms.Path(
		process.hasGenNuMuOrTau
		*process.hasGenNuMuOrTauFilter
		*process.printParticleTree
		)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
		)


process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
		##WR->LNu->LLQQ GEN-SIM files
		'file:/eos/uscms/store/user/skalafut/WR_2000_ToLNu_1000_ToLLQQ_GEN_SIM/WR_M-2000_ToLNu_M-1000_ToLLQQ_13TeV_GEN-SIM/150529_135442/0000/WR_M-2000_ToLNu_M-1000_ToLLQQ_13TeV_GEN_SIM_457.root'

    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)




