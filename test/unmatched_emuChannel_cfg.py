import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEEJJUnmatched")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoEMuChannelSignalUnmatchedModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

from ExoAnalysis.cmsWR.additionalVarParsing_cff import *

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GT, '')


from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")

#################################
#Filters
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.trigFilt = hltHighLevel.clone()
process.trigFilt.HLTPaths = ['HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v*']
process.trigFilt.andOr = True  #if True, then multiple HLT paths will be combined with OR logic

#################################
#Analyzers

## analyze the kinematic distributions of lepton and jet candidates using these modules
process.recoAnalyzerOne = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("recoObjectsAllCuts"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200),
		leptonsOneCollection = cms.InputTag("emuRecoFourObjMass","leptonsOnePassingFourObjMassCut"),#electrons
		leptonsTwoCollection = cms.InputTag("emuRecoFourObjMass","leptonsTwoPassingFourObjMassCut"),#muons
		jetsCollection = cms.InputTag("emuRecoFourObjMass","jetsPassingFourObjMassCut")
		)

process.recoAnalyzerTwo = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("recoObjectsWithPtEtaCuts"),
		leptonsOneCollection = cms.InputTag("emuBareRecoLeptonOne"),#electrons
		leptonsTwoCollection = cms.InputTag("emuBareRecoLeptonTwo"),#muons
		jetsCollection = cms.InputTag("emuBareRecoJetLeptonDrSeparation","jetsPassingDrSeparationCut"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1)
		)

#################################
#Paths
process.unmatchedRecoSignalPath = cms.Path(
		process.trigFilt
		*process.egmGsfElectronIDSequence
		*process.HEEPIDSidebandSequence  #only look for 1 HEEP electron
		*process.wrTunePMuProdSeq
		*process.isHighPtMuSeq
		*process.emuBareRecoParticleSeq
		*process.emuLeadLeptonSeq
		*process.emuDileptonSignalSeq
		*process.emuFilteredRecoDrSeparationSeq
		*process.emuRecoFourObjMassSeq
		*process.recoAnalyzerOne
		)


process.schedule = cms.Schedule(process.unmatchedRecoSignalPath)


process.TFileService = cms.Service("TFileService",
		#fileName = cms.string('/eos/uscms/store/user/skalafut/DoubleEG/analyzed_DoubleEG_skim_signal_region_eejj.root')
		fileName = cms.string(options.output)
	
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)


inputFiles = cms.untracked.vstring(options.files)



process.source = cms.Source( "PoolSource",
	fileNames = inputFiles,
    #inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


