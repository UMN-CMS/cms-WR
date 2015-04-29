import FWCore.ParameterSet.Config as cms


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

## producers to match reco jets and electrons to gen counterparts (gen jets and gen electrons)
matchRecoJetsToGenJetsNoCuts = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedRecoJetsNoCuts"),
		dRforMatching = cms.double(0.3),
		lowLevelCollTag = cms.InputTag("matchGenJetsToGenQuarksNoCuts","matchedGenJetsNoCuts"),
		higherLevelCollTag = cms.InputTag("bareRecoJet")
		)

matchRecoJetsToGenJetsNoCutsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchRecoJetsToGenJetsNoCuts","matchedRecoJetsNoCuts"),
		minNumber = cms.uint32(2)
		)


## match reco electrons to leading and subleading gen electrons from the WR decay chain
## using the next four modules
## no kinematic cuts applied
matchRecoEleToLeadingGenEleNoCuts = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedLeadingRecoEleNoCuts"),
		dRforMatching = cms.double(0.1),
		lowLevelCollTag = cms.InputTag("bareMatchedLeadingGenEle"),
		higherLevelCollTag = cms.InputTag("bareRecoEle")
		)

matchRecoEleToLeadingGenEleNoCutsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchRecoEleToLeadingGenEleNoCuts","matchedLeadingRecoEleNoCuts"),
		minNumber = cms.uint32(1)
		)

matchRecoEleToSubleadingGenEleNoCuts = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedSubleadingRecoEleNoCuts"),
		dRforMatching = cms.double(0.1),
		lowLevelCollTag = cms.InputTag("bareMatchedSubleadingGenEle"),
		higherLevelCollTag = cms.InputTag("bareRecoEle")
		)

matchRecoEleToSubleadingGenEleNoCutsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchRecoEleToSubleadingGenEleNoCuts","matchedSubleadingRecoEleNoCuts"),
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

ptEtaRestrictedMatchedRecoJets = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("matchRecoJetsToGenJetsNoCuts","matchedRecoJetsNoCuts"),
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



