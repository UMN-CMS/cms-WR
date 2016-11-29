import FWCore.ParameterSet.Config as cms

process = cms.Process("AnalyzeWRDecay")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.reformatGenAndRecoCollections_cff')
process.load('ExoAnalysis.cmsWR.genFilterForWrSkims_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('standard') 
options.maxEvents = -1

options.register('channel',
		'',
		VarParsing.VarParsing.multiplicity.singleton,
		VarParsing.VarParsing.varType.string,
		"enable EE or MuMu channel gen filters")

options.parseArguments()

recoLeptons="tempLeptColl"
chargedLeptonId=0
matchedReco=False
hvyNuId=0
if (options.channel=='EE'):
	recoLeptons="slimmedElectrons"
	chargedLeptonId=11
	hvyNuId=9900012
	matchedReco=True
#
if (options.channel=='MuMu'):
	recoLeptons="slimmedMuons"
	chargedLeptonId=13
	hvyNuId=9900014
	matchedReco=False
#

#################################
#Analyzers

## analyze the kinematic distributions of electrons at gen and reco level with these analyzers
## accessing reco leptons matched to GEN leptons currently only works with electrons
process.wrDecayChainAnalyzer = cms.EDAnalyzer('genAndRecoWrAnalyzer',
		treeName = cms.string("genAndMatchedRecoWrDecayNoCuts"),
		leptonPdgId = cms.double(chargedLeptonId),
		minQuarkPdgId = cms.double(1),
		maxQuarkPdgId = cms.double(6),
		firstHeavyParticlePdgId = cms.double(9900024),
		secondHeavyParticlePdgId = cms.double(hvyNuId),
		firstHeavyParticleStatus = cms.double(62),
		saveMatchedRecoInfo = cms.bool(matchedReco),
		cutOnStatusCode = cms.bool(True),
		recoJetCollection = cms.InputTag("slimmedJets"),
		recoLeptonCollection = cms.InputTag(recoLeptons),
		genParticlesCollection = cms.InputTag("genParticlesFormatter"),
		)

#################################
#Paths
process.analyzeWRdecay = cms.Path(
		process.skipGenNuTauSeq
		*process.genAndRecoFormatterSeq
		*process.wrDecayChainAnalyzer
		)

if (options.channel=='EE'):
	process.analyzeWRdecay = cms.Path(process.skipGenNuTauSeq * process.skipGenNuMuSeq * process.genAndRecoFormatterSeq * process.wrDecayChainAnalyzer)

if (options.channel=='MuMu'):
	process.analyzeWRdecay = cms.Path(process.skipGenNuTauSeq * process.skipGenNuEleSeq * process.genAndRecoFormatterSeq * process.wrDecayChainAnalyzer)

process.schedule = cms.Schedule(process.analyzeWRdecay)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string(options.output)
		#fileName = cms.string('analyzedWr.root')
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

inputFiles = cms.untracked.vstring(options.files)

process.source = cms.Source( "PoolSource",
	fileNames = inputFiles,
	#fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/skalafut/public/WR_starting2015/miniAOD/WR_signal_MC/WRToNuEToEEJJ_MW-6000_MNu-3000_TuneCUETP8M1_pythia8_13TeV.root'),
	#inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(options.maxEvents)
    input = cms.untracked.int32(-1)
)


