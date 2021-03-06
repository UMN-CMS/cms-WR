import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEEJJUnmatched")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.recoElectronChannelSignalUnmatchedModules_cff')
process.load('ExoAnalysis.cmsWR.genElectronChannelModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '74X_mcRun2_asymptotic_v2', '')

##jet energy corrections
from ExoAnalysis.cmsWR.JEC_cff import *
JEC_correction(process,'74X_mcRun2_asymptotic_v2')

from ExoAnalysis.cmsWR.heepSelector_cfi import loadHEEPIDSelector
loadHEEPIDSelector(process)
process.load("ExoAnalysis.cmsWR.heepSelector_cfi")

#################################
#Filters
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.trigFilt = hltHighLevel.clone()
process.trigFilt.HLTPaths = ['HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*']
process.trigFilt.andOr = True  #if True, then multiple HLT paths will be combined with OR logic


#################################
#Analyzers

process.genWRAnalyzer = cms.EDAnalyzer('singleParticleAnalyzer',
		treeName = cms.string("genWR"),
		rightHandWsCollection = cms.InputTag("bareMatchedWR")
		)

process.genMatchedParticleAnalyzerAfterFullGenSelection = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("genLeptonsAndJetsWithAllCuts"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200),
		leptonsOneCollection = cms.InputTag("genMatchedFourObjMass","genMatchedLeadingElectronPassingFourObjMassCut"),#lepton with WR mother 
		leptonsTwoCollection = cms.InputTag("genMatchedFourObjMass","genMatchedSubleadingElectronPassingFourObjMassCut"),#lepton with Nu mother
		jetsCollection = cms.InputTag("genMatchedFourObjMass","genJetsPassingFourObjMassCut")
		)


## analyze the kinematic distributions of the reco jets and leptons in WR->ENu->EEJJ evts with these analyzers
process.unmatchedSignalRecoAnalyzerFive = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("signalRecoObjectsWithAllCuts"),
		leptonsCollection = cms.InputTag("recoFourObjMass","leadingLeptonsPassingFourObjMassCut"),
		jetsCollection = cms.InputTag("recoFourObjMass","leadingJetsPassingFourObjMassCut"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.0)
		)

process.unmatchedSignalRecoAnalyzerFivePostHLT = cms.EDAnalyzer('unmatchedAnalyzer',
		treeName = cms.string("signalRecoObjectsWithAllCutsAfterHLT"),
		leptonsCollection = cms.InputTag("recoFourObjMass","leadingLeptonsPassingFourObjMassCut"),
		jetsCollection = cms.InputTag("recoFourObjMass","leadingJetsPassingFourObjMassCut"),
		doDileptonMassCut = cms.bool(True),
		minDileptonMass = cms.double(200.0)
		)

process.genMatchedParticleAnalyzerAfterGenPtEtaDrCuts = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("genLeptonsAndJetsWithPtEtaDrCuts"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsOneCollection = cms.InputTag("genMatchedJetLeptonDrSeparation","matchedLeadingGenElePassingDrSeparationCut"),#lepton with WR mother 
		leptonsTwoCollection = cms.InputTag("genMatchedJetLeptonDrSeparation","matchedSubleadingGenElePassingDrSeparationCut"),#lepton with Nu mother
		jetsCollection = cms.InputTag("genMatchedJetLeptonDrSeparation","genJetsPassingDrSeparationCut")
		)

process.genMatchedParticleAnalyzerAfterGenPtEtaDrCutsPostHLT = cms.EDAnalyzer('unmatchedAnalyzerForMixedLeptonFlavor',
		treeName = cms.string("genLeptonsAndJetsWithPtEtaDrCutsPostHLT"),
		doDileptonMassCut = cms.bool(False),
		minDileptonMass = cms.double(-1),
		leptonsOneCollection = cms.InputTag("genMatchedJetLeptonDrSeparation","matchedLeadingGenElePassingDrSeparationCut"),#lepton with WR mother 
		leptonsTwoCollection = cms.InputTag("genMatchedJetLeptonDrSeparation","matchedSubleadingGenElePassingDrSeparationCut"),#lepton with Nu mother
		jetsCollection = cms.InputTag("genMatchedJetLeptonDrSeparation","genJetsPassingDrSeparationCut")
		)





#################################
#Paths
process.unmatchedRecoSignalPath = cms.Path(
		process.patJetCorrFactorsReapplyJEC
		+process.patJetsReapplyJEC
		+process.trigFilt
		*process.egmGsfElectronIDSequence
		*process.HEEPIDSequence
		*process.bareRecoParticleSeq
		*process.bareRecoDrSeparationSeq
		*process.ptEtaRestrictedSeq
		*process.recoDileptonCandidateSeq
		*process.recoFourObjMassSeq
		*process.unmatchedSignalRecoAnalyzerFive
		)

#use this path to apply the full GEN offline selection, then apply
#the full RECO offline selection without HLT
process.genAndRecoSignalCutsPath = cms.Path(
		#run a quick analyzer to check the number of evts used as input
		process.bareMatchedWRSeq
		*process.genWRAnalyzer
		###identify the gen leptons from the WR decay, and the gen jets
		###matched to the gen quarks from the WR decay
		*process.bareMatchedGenParticleSeq
		*process.bareGenJetSeq
		*process.matchGenJetsToGenQuarksSeq
		##now apply pt, eta, and dR(lepton,jet) cuts to the gen leptons and gen jets
		*process.simultaneousPtEtaCutMatchedObjectsSeq
		#*process.genMatchedParticleAnalyzerTwo
		*process.genMatchedJetLeptonDrSeparationSeq
		#*process.genMatchedParticleAnalyzerTwoPFive	
		##now apply the one lepton pt>60, dilepton mass, and four object mass cuts
		*process.pickGenMatchedEleSeq
		*process.requireGenMatchedHighPtEleSeq
		*process.genMatchedDiLeptonCandidateSeq
		*process.genMatchedFourObjMassSeq
		##run an analyzer to count the number of evts passing the full GEN offline selection
		*process.genMatchedParticleAnalyzerAfterFullGenSelection
		##now apply the full RECO offline selection without HLT and count the number of evts passing
		+process.patJetCorrFactorsReapplyJEC
		+process.patJetsReapplyJEC
		*process.egmGsfElectronIDSequence
		*process.HEEPIDSequence
		*process.bareRecoParticleSeq
		*process.bareRecoDrSeparationSeq
		*process.ptEtaRestrictedSeq
		*process.recoDileptonCandidateSeq
		*process.recoFourObjMassSeq
		*process.unmatchedSignalRecoAnalyzerFive
		)

process.signalCutsHLTEfficiencyPath = cms.Path(
		process.patJetCorrFactorsReapplyJEC
		+process.patJetsReapplyJEC
		+process.egmGsfElectronIDSequence
		*process.HEEPIDSequence
		*process.bareRecoParticleSeq
		*process.bareRecoDrSeparationSeq
		*process.ptEtaRestrictedSeq
		*process.recoDileptonCandidateSeq
		*process.recoFourObjMassSeq
		*process.unmatchedSignalRecoAnalyzerFive
		*process.trigFilt
		*process.unmatchedSignalRecoAnalyzerFivePostHLT
		)

process.genAccCutsHLTEfficiencyPath = cms.Path(
		#identify the gen leptons from the WR decay, and the gen jets
		#matched to the gen quarks from the WR decay
		process.bareMatchedGenParticleSeq
		*process.bareGenJetSeq
		*process.matchGenJetsToGenQuarksSeq
		#now apply pt, eta, and dR(lepton, jet) cuts to the gen leptons and gen jets
		*process.simultaneousPtEtaCutMatchedObjectsSeq
		*process.genMatchedJetLeptonDrSeparationSeq
		*process.genMatchedParticleAnalyzerAfterGenPtEtaDrCuts
		#now apply the trigger filter, and count the number of evts
		#which pass the trigger filter AND the GEN cuts
		*process.trigFilt
		*process.genMatchedParticleAnalyzerAfterGenPtEtaDrCutsPostHLT
		)

process.schedule = cms.Schedule(process.unmatchedRecoSignalPath)

process.TFileService = cms.Service("TFileService",
		fileName = cms.string('analyzed_tree_hltAndAllRecoOfflineCuts_eejjSignalRegion.root')
	
)

process.options = cms.untracked.PSet(
		wantSummary = cms.untracked.bool(True),
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)


#inputFiles = cms.untracked.vstring(options.files)

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
		#'file:/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODBkgndFiles/TTJets_TuneCUETP8M1_13TeV_pythia8_1.root'
		'file:/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/WRToNuEToEEJJ_MW-2000_MNu-1000_TuneCUETP8M1_pythia8_13TeV_5.root'
    ),
    #inputCommands = cms.untracked.vstring(
    #    'keep *'
    #)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


