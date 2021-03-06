import FWCore.ParameterSet.Config as cms
import array

##############################################################
## modules to use on signal samples without GEN matching

bareRecoJet = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedJets"),
		cut = cms.string("")
		)

bareRecoJetFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareRecoJet"),
		minNumber = cms.uint32(2)
		)

bareRecoEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedElectrons"),
		cut = cms.string("")
		)

bareRecoEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareRecoEle"),
		minNumber = cms.uint32(2)
		)

bareRecoParticleSeq = cms.Sequence(bareRecoJet*bareRecoJetFilter*bareRecoEle*bareRecoEleFilter)


## end list of modules to use on signal samples without GEN matching
##############################################################






##############################################################
## modules to use on signal samples with GEN matching

## producers to match reco jets and electrons to gen counterparts (gen jets and gen electrons)
matchRecoJetsToGenJetsNoCuts = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedRecoJetsNoCuts"),
		dRforMatching = cms.double(0.3),
		treeName = cms.string("matchedRecoJetsNoCutsTree"),
		lowLevelCollTag = cms.InputTag("matchGenJetsToGenQuarksNoCuts","matchedGenJetsNoCuts"),
		higherLevelCollTag = cms.InputTag("bareRecoJet")
		)

matchRecoJetsToGenJetsNoCutsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchRecoJetsToGenJetsNoCuts","matchedRecoJetsNoCuts"),
		minNumber = cms.uint32(2)
		)



matchRecoJetsToGenJetsNoCutsNewPath = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedRecoJetsNoCuts"),
		dRforMatching = cms.double(0.3),
		treeName = cms.string("matchedRecoJetsNoCutsTree"),
		lowLevelCollTag = cms.InputTag("matchGenJetsToGenQuarksNoCutsNewPath","matchedGenJetsNoCuts"),
		higherLevelCollTag = cms.InputTag("bareRecoJet")
		)

matchRecoJetsToGenJetsNoCutsNewPathFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchRecoJetsToGenJetsNoCutsNewPath","matchedRecoJetsNoCuts"),
		minNumber = cms.uint32(2)
		)


## producer to match reco jets to gen quarks from the WR decay chain
matchRecoJetsToGenQuarksNoCuts = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedRecoJetsNoCuts"),
		dRforMatching = cms.double(0.35),
		treeName = cms.string("matchedRecoJetsNoCutsTree"),
		lowLevelCollTag = cms.InputTag("bareMatchedGenQuark"),
		higherLevelCollTag = cms.InputTag("bareRecoJet")
		)

matchRecoJetsToGenQuarksNoCutsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchRecoJetsToGenQuarksNoCuts","matchedRecoJetsNoCuts"),
		minNumber = cms.uint32(2)
		)

matchRecoJetsToGenQuarksNoCutsNewPath = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedRecoJetsNoCuts"),
		dRforMatching = cms.double(0.6),
		treeName = cms.string("matchedRecoJetsNoCutsTree"),
		lowLevelCollTag = cms.InputTag("bareMatchedGenQuark"),
		higherLevelCollTag = cms.InputTag("bareRecoJet")
		)

matchRecoJetsToGenQuarksNoCutsNewPathFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchRecoJetsToGenQuarksNoCutsNewPath","matchedRecoJetsNoCuts"),
		minNumber = cms.uint32(2)
		)



## match reco electrons to leading and subleading gen electrons from the WR decay chain
## using the next four modules
## no kinematic cuts applied
matchRecoEleToLeadingGenEleNoCuts = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedLeadingRecoEleNoCuts"),
		dRforMatching = cms.double(0.2),
		treeName = cms.string("matchedLeadingRecoEleNoCutsTree"),
		lowLevelCollTag = cms.InputTag("bareMatchedLeadingGenEle"),
		higherLevelCollTag = cms.InputTag("bareRecoEle")
		)

matchRecoEleToLeadingGenEleNoCutsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchRecoEleToLeadingGenEleNoCuts","matchedLeadingRecoEleNoCuts"),
		minNumber = cms.uint32(1)
		)


matchRecoEleToLeadingGenEleNoCutsNewPath = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedLeadingRecoEleNoCuts"),
		dRforMatching = cms.double(0.2),
		treeName = cms.string("matchedLeadingRecoEleNoCutsTree"),
		lowLevelCollTag = cms.InputTag("bareMatchedLeadingGenEle"),
		higherLevelCollTag = cms.InputTag("bareRecoEle")
		)

matchRecoEleToLeadingGenEleNoCutsNewPathFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchRecoEleToLeadingGenEleNoCutsNewPath","matchedLeadingRecoEleNoCuts"),
		minNumber = cms.uint32(1)
		)


matchRecoEleToSubleadingGenEleNoCuts = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedSubleadingRecoEleNoCuts"),
		dRforMatching = cms.double(0.2),
		treeName = cms.string("matchedSubleadingRecoEleNoCutsTree"),
		lowLevelCollTag = cms.InputTag("bareMatchedSubleadingGenEle"),
		higherLevelCollTag = cms.InputTag("bareRecoEle")
		)

matchRecoEleToSubleadingGenEleNoCutsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchRecoEleToSubleadingGenEleNoCuts","matchedSubleadingRecoEleNoCuts"),
		minNumber = cms.uint32(1)
		)

matchRecoEleToSubleadingGenEleNoCutsNewPath = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedSubleadingRecoEleNoCuts"),
		dRforMatching = cms.double(0.2),
		treeName = cms.string("matchedSubleadingRecoEleNoCutsTree"),
		lowLevelCollTag = cms.InputTag("bareMatchedSubleadingGenEle"),
		higherLevelCollTag = cms.InputTag("bareRecoEle")
		)

matchRecoEleToSubleadingGenEleNoCutsNewPathFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchRecoEleToSubleadingGenEleNoCutsNewPath","matchedSubleadingRecoEleNoCuts"),
		minNumber = cms.uint32(1)
		)


## end reco to gen electron matching without kinematic cuts applied 

matchRecoNoCutsSeq = cms.Sequence(
		matchRecoJetsToGenJetsNoCuts
		*matchRecoJetsToGenJetsNoCutsFilter
		*matchRecoEleToLeadingGenEleNoCuts
		*matchRecoEleToLeadingGenEleNoCutsFilter
		*matchRecoEleToSubleadingGenEleNoCuts
		*matchRecoEleToSubleadingGenEleNoCutsFilter
		)

matchRecoSingleStageJetsNoCutsSeq = cms.Sequence(
		matchRecoJetsToGenQuarksNoCuts
		*matchRecoJetsToGenQuarksNoCutsFilter
		*matchRecoEleToLeadingGenEleNoCuts
		*matchRecoEleToLeadingGenEleNoCutsFilter
		*matchRecoEleToSubleadingGenEleNoCuts
		*matchRecoEleToSubleadingGenEleNoCutsFilter
	
		)

## these next several modules apply pt and eta cuts to the reco objects
## which have been matched to gen objects in the WR decay chain
ptEtaRestrictedMatchedLeadingRecoEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("matchRecoEleToLeadingGenEleNoCuts","matchedLeadingRecoEleNoCuts"),
		cut = cms.string("abs(eta) < 2.5 && pt>60")
		)

ptEtaRestrictedMatchedLeadingRecoEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedMatchedLeadingRecoEle"),
		minNumber = cms.uint32(1)
		)

ptEtaRestrictedMatchedSubleadingRecoEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("matchRecoEleToSubleadingGenEleNoCuts","matchedSubleadingRecoEleNoCuts"),
		cut = cms.string("abs(eta) < 2.5 && pt>40")
		)

ptEtaRestrictedMatchedSubleadingRecoEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedMatchedSubleadingRecoEle"),
		minNumber = cms.uint32(1)
		)

#swap "matchRecoJetsToGenQuarksNoCuts" with "matchedRecoJetsToGenJetsNoCuts" to swap btwn single stage and two stage jet matching
ptEtaRestrictedMatchedRecoJets = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("matchRecoJetsToGenQuarksNoCuts","matchedRecoJetsNoCuts"),
		cut = cms.string("abs(eta) < 2.5 && pt>40")
		)

ptEtaRestrictedMatchedRecoJetsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedMatchedRecoJets"),
		minNumber = cms.uint32(2)
		)

ptEtaRestrictedMatchedRecoSeq = cms.Sequence(
		ptEtaRestrictedMatchedLeadingRecoEle
		*ptEtaRestrictedMatchedLeadingRecoEleFilter
		*ptEtaRestrictedMatchedSubleadingRecoEle
		*ptEtaRestrictedMatchedSubleadingRecoEleFilter
		*ptEtaRestrictedMatchedRecoJets
		*ptEtaRestrictedMatchedRecoJetsFilter
		)
## end modules which apply pt and eta filters to matched reco objects


## these modules apply the dilepton mass cut to matched reco electrons

recoMatchedDiElectronCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("ptEtaRestrictedMatchedLeadingRecoEle ptEtaRestrictedMatchedSubleadingRecoEle"),
		role = cms.string("matchedLeadingEle matchedSubleadingEle"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 200 && daughter(0).pt > daughter(1).pt")
		)

recoMatchedDiElectronCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("recoMatchedDiElectronCandidate"),
		minNumber = cms.uint32(1)
		)

recoMatchedDiElectronCandidateSeq = cms.Sequence(recoMatchedDiElectronCandidate*recoMatchedDiElectronCandidateFilter)

## end modules which apply dilepton mass cut to reco electrons matched to gen electrons


## end list of modules to use on signal samples with GEN matching
##############################################################


