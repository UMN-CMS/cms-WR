import FWCore.ParameterSet.Config as cms

##filters to select gen particles matched to WR decay products via pdgId
##heavyNu muon flavor has pdgId = 9900014
##WR has pdgId = 9900024

#this filter looks for leptons whose mother is a WR, has |eta| < 2.4, and pt>60
mumuChnlPtEtaFilteredMatchedLeadingGenLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("genParticles"),
		cut = cms.string("abs(pdgId) == 13 && abs(mother(0).pdgId) == 9900024 && abs(eta) < 2.4 && pt>60")
		)

mumuChnlPtEtaFilteredMatchedLeadingGenLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("mumuChnlPtEtaFilteredMatchedLeadingGenLepton"),
		minNumber = cms.uint32(1)
		)

#this filter looks for muons whose mother is a heavy Nu, has |eta| < 2.4, and pt>40
mumuChnlPtEtaFilteredMatchedSubleadingGenLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("genParticles"),
		cut = cms.string("abs(pdgId) == 13 && abs(mother(0).pdgId) == 9900014 && abs(eta) < 2.4 && pt>40")
		)

mumuChnlPtEtaFilteredMatchedSubleadingGenLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("mumuChnlPtEtaFilteredMatchedSubleadingGenLepton"),
		minNumber = cms.uint32(1)
		)

#this filter looks for quarks whose real mother is a heavy Nu (virtuals are not tracked in Pythia), |eta| < 2.5, and pt>40
mumuChnlPtEtaFilteredMatchedGenQuark = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("genParticles"),
		cut = cms.string("abs(pdgId) < 7 && abs(mother(0).pdgId) == 9900014 && abs(eta) < 2.5 && pt>40")
		)

mumuChnlPtEtaFilteredMatchedGenQuarkFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("mumuChnlPtEtaFilteredMatchedGenQuark"),
		minNumber = cms.uint32(2)
		)

mumuChnlPtEtaFilteredMatchedGenParticleSeq = cms.Sequence(
		mumuChnlPtEtaFilteredMatchedLeadingGenLepton
		*mumuChnlPtEtaFilteredMatchedLeadingGenLeptonFilter
		*mumuChnlPtEtaFilteredMatchedSubleadingGenLepton
		*mumuChnlPtEtaFilteredMatchedSubleadingGenLeptonFilter
		*mumuChnlPtEtaFilteredMatchedGenQuark
		*mumuChnlPtEtaFilteredMatchedGenQuarkFilter
		)


##filters on dilepton mass, and the producer necessary to facilitate this filter
mumuChnlMatchedDiLeptonCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("mumuChnlPtEtaFilteredMatchedLeadingGenLepton mumuChnlPtEtaFilteredMatchedSubleadingGenLepton"),
		role = cms.string("matchedLeadingLepton matchedSubleadingLepton"),
		checkCharge = cms.bool(False),
		#cut = cms.string("mass > 200 && daughter(0).pt > daughter(1).pt")
		#when the Nu/WR mass ratio is very small the lepton with Nu mother
		#will generally be harder than the lepton with WR mother
		cut = cms.string("mass > 200")
		)

mumuChnlMatchedDiLeptonCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("mumuChnlMatchedDiLeptonCandidate"),
		minNumber = cms.uint32(1)
		)

mumuChnlMatchedDiLeptonCandidateSeq = cms.Sequence(mumuChnlMatchedDiLeptonCandidate*mumuChnlMatchedDiLeptonCandidateFilter)

##skip evts where two quarks are not found with |eta| < 2.5, pt>40, and are separated from both leptons by some dR value 
##in addition, the leptons must be separated from each other by some other dR value 
mumuChnlMatchedQuarkLeptonDrSeparation = cms.EDProducer("applyLeptonJetDrCutMixedLeptonFlavor",
		inputJetsCollTag = cms.InputTag("mumuChnlPtEtaFilteredMatchedGenQuark"),
		inputLeptonsOneCollTag = cms.InputTag("mumuChnlPtEtaFilteredMatchedLeadingGenLepton"),
		inputLeptonsTwoCollTag = cms.InputTag("mumuChnlPtEtaFilteredMatchedSubleadingGenLepton"),
		minDileptonMassCut = cms.double(200),
		minDrSeparation = cms.double(0.4),#lepton jet separation
		outputJetsCollectionName = cms.string("drSeparatedJetsPassingPtEtaFilters"),
		minLeptonDrSeparation = cms.double(0.4) #the two selected leptons must be separated by at least this much in dR
		)

mumuChnlMatchedQuarkLeptonDrSeparationFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("mumuChnlMatchedQuarkLeptonDrSeparation","drSeparatedJetsPassingPtEtaFilters"),
		minNumber = cms.uint32(2)
		)

mumuChnlMatchedQuarkLeptonDrSepSeq = cms.Sequence(mumuChnlMatchedQuarkLeptonDrSeparation*mumuChnlMatchedQuarkLeptonDrSeparationFilter)


##apply four object mass cut
##if multiple four object candidates are found in one event, then the constituents of each four object will be saved to the output collections
mumuChnlMatchedFourObjectCandidate = cms.EDProducer("applyFourObjMassCutTwoInputLeptonColls",
		outputJetsCollectionName = cms.string("matchedJetsPassingAllCuts"),
		inputJetsCollTag = cms.InputTag("mumuChnlMatchedQuarkLeptonDrSeparation","drSeparatedJetsPassingPtEtaFilters"),
		outputLeptonsOneCollectionName = cms.string("matchedLeadingGenLeptonPassingAllCuts"),
		inputLeptonsOneCollTag = cms.InputTag("mumuChnlPtEtaFilteredMatchedLeadingGenLepton"),
		outputLeptonsTwoCollectionName = cms.string("matchedSubleadingGenLeptonPassingAllCuts"),
		inputLeptonsTwoCollTag = cms.InputTag("mumuChnlPtEtaFilteredMatchedSubleadingGenLepton"),
		minFourObjMassCut = cms.double(600),
		minDileptonMassCut = cms.double(200),
		minLeptonDrSeparation = cms.double(0.4)
		)

mumuChnlMatchedFourObjectCandidateLeptonOneCheck = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("mumuChnlMatchedFourObjectCandidate","matchedLeadingGenLeptonPassingAllCuts"),
		minNumber = cms.uint32(1)
		)

mumuChnlMatchedFourObjectCandidateLeptonTwoCheck = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("mumuChnlMatchedFourObjectCandidate","matchedSubleadingGenLeptonPassingAllCuts"),
		minNumber = cms.uint32(1)
		)

mumuChnlMatchedFourObjectCandidateJetsCheck = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("mumuChnlMatchedFourObjectCandidate","matchedJetsPassingAllCuts"),
		minNumber = cms.uint32(1)
		)

mumuChnlMatchedFourObjectCandidateSeq = cms.Sequence(
		mumuChnlMatchedFourObjectCandidate
		*mumuChnlMatchedFourObjectCandidateLeptonOneCheck
		*mumuChnlMatchedFourObjectCandidateLeptonTwoCheck
		*mumuChnlMatchedFourObjectCandidateJetsCheck
		)

