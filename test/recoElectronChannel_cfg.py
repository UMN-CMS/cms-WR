import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOEEJJ")

## load the filters, producers, and sequences defined in other config file fragments
process.load('ExoAnalysis.cmsWR.genElectronChannelModules_cff')
process.load('ExoAnalysis.cmsWR.recoElectronChannelModules_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('standard') 
options.maxEvents = -1
options.parseArguments()

#################################
#Analyzers

## analyze the kinematic distributions of the reco jets and leptons in WR->ENu->EEJJ evts with these analyzers 

#for swapping btwn two stage and single stage jet matching:
##change "matchGenJetsToGenQuarksNoCuts","matchedGenJetsNoCuts" to "bareMatchedGenQuark"
##and "matchRecoJetsToGenJetsNoCuts" to "matchRecoJetsToGenQuarksNoCuts"

process.matchedRecoAnalyzerOne = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedRecoObjectsNoCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		saveGenMatched = cms.bool(True),
		leadingLeptonCollection = cms.InputTag("matchRecoEleToLeadingGenEleNoCuts","matchedLeadingRecoEleNoCuts"),
		subleadingLeptonCollection = cms.InputTag("matchRecoEleToSubleadingGenEleNoCuts","matchedSubleadingRecoEleNoCuts"),
		quarkCollection = cms.InputTag("matchRecoJetsToGenQuarksNoCuts","matchedRecoJetsNoCuts"),
		matchingLeadingLeptonCollection = cms.InputTag("bareMatchedLeadingGenEle"),
		matchingSubleadingLeptonCollection = cms.InputTag("bareMatchedSubleadingGenEle"),
		matchingQuarkCollection = cms.InputTag("bareMatchedGenQuark")
		)


process.matchedRecoAnalyzerTwo = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedRecoObjectsWithPtEtaCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		saveGenMatched = cms.bool(True),
		leadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingRecoEle"),
		subleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingRecoEle"),
		quarkCollection = cms.InputTag("ptEtaRestrictedMatchedRecoJets"),
		matchingLeadingLeptonCollection = cms.InputTag("bareMatchedLeadingGenEle"),
		matchingSubleadingLeptonCollection = cms.InputTag("bareMatchedSubleadingGenEle"),
		matchingQuarkCollection = cms.InputTag("bareMatchedGenQuark")
		)

process.matchedRecoAnalyzerThree = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedRecoObjectsWithPtEtaAndDileptonMassCuts"),
		doDeltaRcut = cms.bool(False),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(-1),
		saveGenMatched = cms.bool(True),
		leadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingRecoEle"),
		subleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingRecoEle"),
		quarkCollection = cms.InputTag("ptEtaRestrictedMatchedRecoJets"),
		matchingLeadingLeptonCollection = cms.InputTag("bareMatchedLeadingGenEle"),
		matchingSubleadingLeptonCollection = cms.InputTag("bareMatchedSubleadingGenEle"),
		matchingQuarkCollection = cms.InputTag("bareMatchedGenQuark")
		)

process.matchedRecoAnalyzerFour = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedRecoObjectsWithPtEtaDileptonMassAndDrCuts"),
		doDeltaRcut = cms.bool(True),
		doFourObjMassCut = cms.bool(False),
		minFourObjMass = cms.double(-1),
		minDeltaRforLeptonJetExclusion = cms.double(0.4),
		saveGenMatched = cms.bool(True),
		leadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingRecoEle"),
		subleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingRecoEle"),
		quarkCollection = cms.InputTag("ptEtaRestrictedMatchedRecoJets"),
		matchingLeadingLeptonCollection = cms.InputTag("bareMatchedLeadingGenEle"),
		matchingSubleadingLeptonCollection = cms.InputTag("bareMatchedSubleadingGenEle"),
		matchingQuarkCollection = cms.InputTag("bareMatchedGenQuark")
		)


process.matchedRecoAnalyzerFive = cms.EDAnalyzer('matchedAnalyzer',
		treeName = cms.string("matchedRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts"),
		doDeltaRcut = cms.bool(True),
		doFourObjMassCut = cms.bool(True),
		minFourObjMass = cms.double(600),
		minDeltaRforLeptonJetExclusion = cms.double(0.4),
		saveGenMatched = cms.bool(True),
		leadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedLeadingRecoEle"),
		subleadingLeptonCollection = cms.InputTag("ptEtaRestrictedMatchedSubleadingRecoEle"),
		quarkCollection = cms.InputTag("ptEtaRestrictedMatchedRecoJets"),
		matchingLeadingLeptonCollection = cms.InputTag("bareMatchedLeadingGenEle"),
		matchingSubleadingLeptonCollection = cms.InputTag("bareMatchedSubleadingGenEle"),
		matchingQuarkCollection = cms.InputTag("bareMatchedGenQuark")
		)


#################################
#Paths
process.matchedRecoTwoStageJetsPath = cms.Path(
		process.bareGenJet
		* process.bareGenJetFilter
		* process.bareMatchedGenParticleSeq
		* process.matchGenJetsToGenQuarksNoCuts
		* process.matchGenJetsToGenQuarksNoCutsFilter
		* process.bareRecoParticleSeq
		* process.matchRecoNoCutsSeq
		* process.matchedRecoAnalyzerOne
		* process.ptEtaRestrictedMatchedRecoSeq
		* process.matchedRecoAnalyzerTwo
		* process.recoMatchedDiElectronCandidateSeq
		* process.matchedRecoAnalyzerThree
		* process.matchedRecoAnalyzerFour
		* process.matchedRecoAnalyzerFive

		)

process.matchedRecoSingleStageJetsPath = cms.Path(
		process.bareMatchedGenParticleSeq
		* process.bareRecoParticleSeq
		* process.matchRecoSingleStageJetsNoCutsSeq
		* process.matchedRecoAnalyzerOne
		* process.ptEtaRestrictedMatchedRecoSeq
		* process.matchedRecoAnalyzerTwo
		* process.recoMatchedDiElectronCandidateSeq
		* process.matchedRecoAnalyzerThree
		* process.matchedRecoAnalyzerFour
		* process.matchedRecoAnalyzerFive
		)


process.onlyMatchTwoStageJetsPath = cms.Path(
		process.bareGenJet
		* process.bareGenJetFilter
		* process.bareMatchedGenParticleSeq
		* process.matchGenJetsToGenQuarksNoCutsNewPath
		* process.matchGenJetsToGenQuarksNoCutsNewPathFilter
		* process.bareRecoJet
		* process.bareRecoJetFilter
		* process.matchRecoJetsToGenJetsNoCutsNewPath
		* process.matchRecoJetsToGenJetsNoCutsNewPathFilter
		)

process.onlyMatchSingleStageJetsPath = cms.Path(
		process.bareMatchedGenParticleSeq
		* process.bareRecoJet
		* process.bareRecoJetFilter
		* process.matchRecoJetsToGenQuarksNoCutsNewPath
		* process.matchRecoJetsToGenQuarksNoCutsNewPathFilter
		)

process.onlyMatchLeadingElePath = cms.Path(
		process.bareMatchedGenParticleSeq
		* process.bareRecoEle
		* process.bareRecoEleFilter
		* process.matchRecoEleToLeadingGenEleNoCutsNewPath
		* process.matchRecoEleToLeadingGenEleNoCutsNewPathFilter
		)

process.onlyMatchSubleadingElePath = cms.Path(
		process.bareMatchedGenParticleSeq
		* process.bareRecoEle
		* process.bareRecoEleFilter
		* process.matchRecoEleToSubleadingGenEleNoCutsNewPath
		* process.matchRecoEleToSubleadingGenEleNoCutsNewPathFilter
		)

process.schedule = cms.Schedule(process.matchedRecoSingleStageJetsPath)


process.TFileService = cms.Service("TFileService",
		#fileName = cms.string('analysis_recoElectronChannel_two_stage_matching_for_jets.root')
		#fileName = cms.string('analysis_recoElectronChannel_single_stage_matching_for_jets_centrally_produced_signal_MWr_2600_MNu_1300.root')
		#fileName = cms.string('analysis_ttbar_bkgndElectronChannel.root')
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


