import FWCore.ParameterSet.Config as cms

printParticleTree = cms.EDAnalyzer("ParticleListDrawer",
		maxEventsToPrint = cms.untracked.int32(-1),
		printVertex = cms.untracked.bool(False),
		printOnlyHardInteraction = cms.untracked.bool(False),
		src = cms.InputTag("genParticles")
		#src = cms.InputTag("prunedGenParticles")
		)

hasGenNuMuOrTau = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("genParticles"),
		cut = cms.string("abs(pdgId) > 9900012 && abs(pdgId) < 9900024")
		)

hasGenNuMuOrTauFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("hasGenNuMuOrTau"),
		minNumber = cms.uint32(1)
		)

bareGenJet = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedGenJets"),
		#src = cms.InputTag("ak4GenJets"),
		cut = cms.string("")
		)

bareGenJetFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareGenJet"),
		minNumber = cms.uint32(2)
		)

bareGenEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		#src = cms.InputTag("genParticles"),
		cut = cms.string("abs(pdgId) == 13")
		)

bareGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareGenEle"),
		minNumber = cms.uint32(2)
		)

bareGenParticleSeq = cms.Sequence(bareGenJet*bareGenJetFilter*bareGenEle*bareGenEleFilter)

## producers to find GenJets matched to Gen quarks

matchGenJetsToGenQuarksNoCuts = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedGenJetsNoCuts"),
		dRforMatching = cms.double(0.6),
		treeName = cms.string("matchedGenJetsNoCutsTree"),
		lowLevelCollTag = cms.InputTag("bareMatchedGenQuark"),
		higherLevelCollTag = cms.InputTag("bareGenJet")
		)

matchGenJetsToGenQuarksNoCutsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchGenJetsToGenQuarksNoCuts","matchedGenJetsNoCuts"),
		minNumber = cms.uint32(2)
		)

matchGenJetsToGenQuarksNoCutsNewPath = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedGenJetsNoCuts"),
		dRforMatching = cms.double(0.6),
		treeName = cms.string("matchedGenJetsNoCutsTree"),
		lowLevelCollTag = cms.InputTag("bareMatchedGenQuark"),
		higherLevelCollTag = cms.InputTag("bareGenJet")
		)

matchGenJetsToGenQuarksNoCutsNewPathFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchGenJetsToGenQuarksNoCutsNewPath","matchedGenJetsNoCuts"),
		minNumber = cms.uint32(2)
		)




##filters to select gen particles matched to WR decay products via pdgId
##heavyNu pdgId = 9900014, WR pdgId = 9900024

#this filter looks for electrons whose mother is a WR
bareMatchedLeadingGenEle = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("genParticles"),
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) == 13 && abs(mother(0).pdgId) == 9900024")
		)

bareMatchedLeadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareMatchedLeadingGenEle"),
		minNumber = cms.uint32(1)
		)

#this filter looks for electrons whose mother is a heavy Nu
bareMatchedSubleadingGenEle = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("genParticles"),
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) == 13 && abs(mother(0).pdgId) == 9900014")
		)

bareMatchedSubleadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareMatchedSubleadingGenEle"),
		minNumber = cms.uint32(1)
		)

#this filter looks for quarks whose real mother is a heavy Nu (virtuals are not tracked in Pythia) 
bareMatchedGenQuark = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("genParticles"),
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) < 7 && abs(mother(0).pdgId) == 9900014")
		)

bareMatchedGenQuarkFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareMatchedGenQuark"),
		minNumber = cms.uint32(2)
		)

bareMatchedGenParticleSeq = cms.Sequence(bareMatchedLeadingGenEle*bareMatchedLeadingGenEleFilter*bareMatchedSubleadingGenEle*
		bareMatchedSubleadingGenEleFilter*bareMatchedGenQuark*bareMatchedGenQuarkFilter)


##filters to cut out gen particles which could never hit the tracker based on eta trajectory
##use these modules to determine efficiency of signal skim (skim applies pT cuts to leptons and possibly jets)
##skim efficiency on signal = # evts passing reco cuts / # evts before reco cuts with gen objects falling within tracker coverage
etaRestrictedMatchedGenQuark = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareMatchedGenQuark"),
		cut = cms.string("abs(eta) < 2.5")
		)
etaRestrictedMatchedGenQuarkFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("etaRestrictedMatchedGenQuark"),
		minNumber = cms.uint32(2)
		)

etaRestrictedMatchedGenLeadingLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareMatchedLeadingGenEle"),
		cut = cms.string("abs(eta) < 2.4")
		)
etaRestrictedMatchedGenLeadingLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("etaRestrictedMatchedGenLeadingLepton"),
		minNumber = cms.uint32(1)
		)

etaRestrictedMatchedGenSubleadingLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareMatchedSubleadingGenEle"),
		cut = cms.string("abs(eta) < 2.4")
		)
etaRestrictedMatchedGenSubleadingLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("etaRestrictedMatchedGenSubleadingLepton"),
		minNumber = cms.uint32(1)
		)

etaRestrictedMatchedGenParticleSeq = cms.Sequence(
		etaRestrictedMatchedGenQuark
		*etaRestrictedMatchedGenQuarkFilter
		*etaRestrictedMatchedGenLeadingLepton
		*etaRestrictedMatchedGenLeadingLeptonFilter
		*etaRestrictedMatchedGenSubleadingLepton
		*etaRestrictedMatchedGenSubleadingLeptonFilter
		)


##filters on pt and eta of gen leptons and jets
ptEtaRestrictedGenJet = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("slimmedGenJets"),
		src = cms.InputTag("bareGenJet"),
		cut = cms.string("pt>40 && abs(eta) < 2.5")
		)

ptEtaRestrictedGenJetFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedGenJet"),
		minNumber = cms.uint32(2)
		)

ptEtaRestrictedSubleadingGenEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareGenEle"),
		cut = cms.string("pt>40 && abs(eta) < 2.4")
		)

ptEtaRestrictedSubleadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedSubleadingGenEle"),
		minNumber = cms.uint32(2)
		)

ptEtaRestrictedLeadingGenEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("ptEtaRestrictedSubleadingGenEle"),
		cut = cms.string("pt>60")
		)

ptEtaRestrictedLeadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedLeadingGenEle"),
		minNumber = cms.uint32(1)
		)

ptEtaConstrainedSeq = cms.Sequence(ptEtaRestrictedGenJet*ptEtaRestrictedGenJetFilter*ptEtaRestrictedSubleadingGenEle
		*ptEtaRestrictedSubleadingGenEleFilter*ptEtaRestrictedLeadingGenEle*ptEtaRestrictedLeadingGenEleFilter)


#these filters and sequence look for gen quarks and leptons which are matched (by pdgId) to particles in the WR
#decay chain

ptEtaRestrictedMatchedGenQuark = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("etaRestrictedMatchedGenQuark"),
		#src = cms.InputTag("bareMatchedGenQuark"),
		cut = cms.string("pt>40 && abs(eta) < 2.5")
		)

ptEtaRestrictedMatchedGenQuarkFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedMatchedGenQuark"),
		minNumber = cms.uint32(2)
		)

ptEtaRestrictedMatchedSubleadingGenEle = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("bareMatchedSubleadingGenEle"),
		src = cms.InputTag("etaRestrictedMatchedGenSubleadingLepton"),
		cut = cms.string("pt>40 && abs(eta) < 2.4")
		)

ptEtaRestrictedMatchedSubleadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedMatchedSubleadingGenEle"),
		minNumber = cms.uint32(1)
		)

ptEtaRestrictedMatchedLeadingGenEle = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("bareMatchedLeadingGenEle"),
		src = cms.InputTag("etaRestrictedMatchedGenLeadingLepton"),
		cut = cms.string("pt>60 && abs(eta) < 2.4")
		)

ptEtaRestrictedMatchedLeadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedMatchedLeadingGenEle"),
		minNumber = cms.uint32(1)
		)

matchedPtEtaConstrainedSeq = cms.Sequence(ptEtaRestrictedMatchedGenQuark*ptEtaRestrictedMatchedGenQuarkFilter
		*ptEtaRestrictedMatchedSubleadingGenEle*ptEtaRestrictedMatchedSubleadingGenEleFilter
		*ptEtaRestrictedMatchedLeadingGenEle*ptEtaRestrictedMatchedLeadingGenEleFilter)




##filters on dilepton mass, and the producer necessary to facilitate this filter
genDiElectronCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("ptEtaRestrictedLeadingGenEle ptEtaRestrictedSubleadingGenEle"),
		role = cms.string("leading subleading"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 200 && daughter(0).pt > daughter(1).pt")
		)

genDiElectronCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genDiElectronCandidate"),
		minNumber = cms.uint32(1)
		)

genDiElectronCandidateSeq = cms.Sequence(genDiElectronCandidate*genDiElectronCandidateFilter)

genMatchedDiElectronCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("ptEtaRestrictedMatchedLeadingGenEle ptEtaRestrictedMatchedSubleadingGenEle"),
		role = cms.string("matchedLeadingEle matchedSubleadingEle"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 200 && daughter(0).pt > daughter(1).pt")
		)

genMatchedDiElectronCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genMatchedDiElectronCandidate"),
		minNumber = cms.uint32(1)
		)

genMatchedDiElectronCandidateSeq = cms.Sequence(genMatchedDiElectronCandidate*genMatchedDiElectronCandidateFilter)



