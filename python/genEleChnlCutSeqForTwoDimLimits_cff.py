import FWCore.ParameterSet.Config as cms

##filters to select gen particles matched to WR decay products via pdgId
##heavyNu electron flavor has pdgId = 9900012
##WR has pdgId = 9900024

#this filter looks for leptons whose mother is a WR, has |eta| < 2.5, and pt>60
eeChnlPtEtaFilteredMatchedLeadingGenLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("genParticles"),
		cut = cms.string("abs(pdgId) == 11 && abs(mother(0).pdgId) == 9900024 && abs(eta) < 2.5 && pt>60")
		)

eeChnlPtEtaFilteredMatchedLeadingGenLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("eeChnlPtEtaFilteredMatchedLeadingGenLepton"),
		minNumber = cms.uint32(1)
		)

#this filter looks for electrons whose mother is a heavy Nu, has |eta| < 2.5, and pt>40
eeChnlPtEtaFilteredMatchedSubleadingGenLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("genParticles"),
		cut = cms.string("abs(pdgId) == 11 && abs(mother(0).pdgId) == 9900012 && abs(eta) < 2.5 && pt>40")
		)

eeChnlPtEtaFilteredMatchedSubleadingGenLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("eeChnlPtEtaFilteredMatchedSubleadingGenLepton"),
		minNumber = cms.uint32(1)
		)

#this filter looks for quarks whose real mother is a heavy Nu (virtuals are not tracked in Pythia), |eta| < 2.5, and pt>40
eeChnlPtEtaFilteredMatchedGenQuark = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("genParticles"),
		cut = cms.string("abs(pdgId) < 7 && abs(mother(0).pdgId) == 9900012 && abs(eta) < 2.5 && pt>40")
		)

eeChnlPtEtaFilteredMatchedGenQuarkFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("eeChnlPtEtaFilteredMatchedGenQuark"),
		minNumber = cms.uint32(2)
		)

eeChnlPtEtaFilteredMatchedGenParticleSeq = cms.Sequence(
		eeChnlPtEtaFilteredMatchedLeadingGenLepton
		*eeChnlPtEtaFilteredMatchedLeadingGenLeptonFilter
		*eeChnlPtEtaFilteredMatchedSubleadingGenLepton
		*eeChnlPtEtaFilteredMatchedSubleadingGenLeptonFilter
		*eeChnlPtEtaFilteredMatchedGenQuark
		*eeChnlPtEtaFilteredMatchedGenQuarkFilter
		)


##filters on dilepton mass, and the producer necessary to facilitate this filter
eeChnlMatchedDiLeptonCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("eeChnlPtEtaFilteredMatchedLeadingGenLepton eeChnlPtEtaFilteredMatchedSubleadingGenLepton"),
		role = cms.string("matchedLeadingLepton matchedSubleadingLepton"),
		checkCharge = cms.bool(False),
		#cut = cms.string("mass > 200 && daughter(0).pt > daughter(1).pt")
		#when the Nu/WR mass ratio is very small the lepton with Nu mother
		#will generally be harder than the lepton with WR mother
		cut = cms.string("mass > 200")
		)

eeChnlMatchedDiLeptonCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("eeChnlMatchedDiLeptonCandidate"),
		minNumber = cms.uint32(1)
		)

eeChnlMatchedDiLeptonCandidateSeq = cms.Sequence(eeChnlMatchedDiLeptonCandidate*eeChnlMatchedDiLeptonCandidateFilter)

##skip evts where two quarks are not found with |eta| < 2.5, pt>40, and are separated from both leptons by some dR value 
##in addition, the leptons must be separated from each other by some other dR value 
eeChnlMatchedQuarkLeptonDrSeparation = cms.EDProducer("applyLeptonJetDrCutMixedLeptonFlavor",
		inputJetsCollTag = cms.InputTag("eeChnlPtEtaFilteredMatchedGenQuark"),
		inputLeptonsOneCollTag = cms.InputTag("eeChnlPtEtaFilteredMatchedLeadingGenLepton"),
		inputLeptonsTwoCollTag = cms.InputTag("eeChnlPtEtaFilteredMatchedSubleadingGenLepton"),
		minDileptonMassCut = cms.double(200),
		minDrSeparation = cms.double(0.4),#lepton jet separation
		outputJetsCollectionName = cms.string("drSeparatedJetsPassingPtEtaFilters"),
		minLeptonDrSeparation = cms.double(0.4) #the two selected leptons must be separated by at least this much in dR
		)

eeChnlMatchedQuarkLeptonDrSeparationFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("eeChnlMatchedQuarkLeptonDrSeparation","drSeparatedJetsPassingPtEtaFilters"),
		minNumber = cms.uint32(2)
		)

eeChnlMatchedQuarkLeptonDrSepSeq = cms.Sequence(eeChnlMatchedQuarkLeptonDrSeparation*eeChnlMatchedQuarkLeptonDrSeparationFilter)


##apply four object mass cut
##if multiple four object candidates are found in one event, then the constituents of each four object will be saved to the output collections
eeChnlMatchedFourObjectCandidate = cms.EDProducer("applyFourObjMassCutTwoInputLeptonColls",
		outputJetsCollectionName = cms.string("matchedJetsPassingAllCuts"),
		inputJetsCollTag = cms.InputTag("eeChnlMatchedQuarkLeptonDrSeparation","drSeparatedJetsPassingPtEtaFilters"),
		outputLeptonsOneCollectionName = cms.string("matchedLeadingGenLeptonPassingAllCuts"),
		inputLeptonsOneCollTag = cms.InputTag("eeChnlPtEtaFilteredMatchedLeadingGenLepton"),
		outputLeptonsTwoCollectionName = cms.string("matchedSubleadingGenLeptonPassingAllCuts"),
		inputLeptonsTwoCollTag = cms.InputTag("eeChnlPtEtaFilteredMatchedSubleadingGenLepton"),
		minFourObjMassCut = cms.double(600),
		minDileptonMassCut = cms.double(200),
		minLeptonDrSeparation = cms.double(0.4)
		)

eeChnlMatchedFourObjectCandidateLeptonOneCheck = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("eeChnlMatchedFourObjectCandidate","matchedLeadingGenLeptonPassingAllCuts"),
		minNumber = cms.uint32(1)
		)

eeChnlMatchedFourObjectCandidateLeptonTwoCheck = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("eeChnlMatchedFourObjectCandidate","matchedSubleadingGenLeptonPassingAllCuts"),
		minNumber = cms.uint32(1)
		)

eeChnlMatchedFourObjectCandidateJetsCheck = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("eeChnlMatchedFourObjectCandidate","matchedJetsPassingAllCuts"),
		minNumber = cms.uint32(1)
		)

eeChnlMatchedFourObjectCandidateSeq = cms.Sequence(
		eeChnlMatchedFourObjectCandidate
		*eeChnlMatchedFourObjectCandidateLeptonOneCheck
		*eeChnlMatchedFourObjectCandidateLeptonTwoCheck
		*eeChnlMatchedFourObjectCandidateJetsCheck
		)

