import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEMu")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoEMuChannelSidebandUnmatchedModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

from ExoAnalysis.cmsWR.additionalVarParsing_cff import *

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GT, '')

##jet energy corrections
from ExoAnalysis.cmsWR.JEC_cff import *
JEC_correction(process, options.GT)

# import the HEEP selection modules and sequences
from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")


#################################
#Filters
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.trigFilt = hltHighLevel.clone()
#process.trigFilt.HLTPaths = ['HLT_Mu45_eta2p1_v*','HLT_Mu50_v*','HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*']
process.trigFilt.HLTPaths = ['HLT_Mu45_eta2p1_v*']
#process.trigFilt.HLTPaths = ['HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v*']
process.trigFilt.andOr = True  #if True, then multiple HLT paths will be combined with OR logic

#################################
#Analyzers

## analyze the kinematic distributions of lepton candidates using these modules
process.recoAnalyzerOne = cms.EDAnalyzer('emuAnalyzer',
		treeName = cms.string("recoObjectsNoCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		minDr = cms.double(0.4),
		leptonsOneCollection = cms.InputTag("electronCheck"),#electrons
		leptonsTwoCollection = cms.InputTag("muonCheck")#muons
		)

#################################
#Paths
process.unmatchedBkgndRecoPath = cms.Path(
		process.patJetCorrFactorsReapplyJEC
		+process.patJetsReapplyJEC
		#process.trigFilt
		+process.egmGsfElectronIDSequence
		*process.HEEPIDSidebandSequence  #only look for 1 HEEP electron
		*process.wrTunePMuProdSeq
		*process.isHighPtMuSeq
		*process.checkEMuSeq
		*process.recoAnalyzerOne
		)
process.schedule = cms.Schedule(process.unmatchedBkgndRecoPath)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string(options.output)

)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True)
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


