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
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

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
		process.trigFilt
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
		process.egmGsfElectronIDSequence
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

#process.schedule = cms.Schedule(process.unmatchedRecoSignalPath)
process.schedule = cms.Schedule(process.signalCutsHLTEfficiencyPath,process.genAccCutsHLTEfficiencyPath)


process.TFileService = cms.Service("TFileService",
		fileName = cms.string('analyzed_tree_checkHLTEfficiencyRecoFullOfflineAndGenPtEtaDr_eejjSignalRegion.root')
	
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


