import FWCore.ParameterSet.Config as cms

process = cms.Process("CheckWRDecay")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.genElectronChannelModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

#import FWCore.ParameterSet.VarParsing as VarParsing
#options = VarParsing.VarParsing('standard') 
#options.maxEvents = -1
#options.parseArguments()


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
		process.bareMatchedWRSeq
		*process.genWRAnalyzerOne
		)
process.schedule = cms.Schedule(process.checkWRdecay)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('file:genWrKinematics.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)


process.source = cms.Source( "PoolSource",
	fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuEToEEJJ_MW-1400_MNu-700-TuneCUETP8M1_pythia8_13TeV_1.root'),
	#inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


