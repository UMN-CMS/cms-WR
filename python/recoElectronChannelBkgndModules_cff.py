import FWCore.ParameterSet.Config as cms

## transform the objects in slimmedJets and slimmedElectrons into reco::Candidate objects
bareRecoJet = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedJets"),
		cut = cms.string("")
		)

bareRecoJetFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareRecoJet"),
		minNumber = cms.uint32(2)
		)

bareRecoLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedElectrons"),
		cut = cms.string("")
		)

bareRecoLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareRecoLepton"),
		minNumber = cms.uint32(2)
		)

bareRecoParticleSeq = cms.Sequence(bareRecoJet*bareRecoJetFilter*bareRecoLepton*bareRecoLeptonFilter)

## apply pT and eta cuts to the bare reco::Candidate leptons and jets
ptEtaRestrictedRecoLeptons = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareRecoLepton"),
		cut = cms.string("abs(eta) < 2.5 && pt>40")
		)
ptEtaRestrictedRecoLeptonsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		minNumber = cms.uint32(2)
		)

ptEtaRestrictedLeadRecoLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		cut = cms.string("pt>60")
		)
ptEtaRestrictedLeadRecoLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedLeadRecoLepton"),
		minNumber = cms.uint32(1)
		)

ptEtaRestrictedRecoJets = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareRecoJet"),
		cut = cms.string("abs(eta) < 2.5 && pt>40")
		)
ptEtaRestrictedRecoJetsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedRecoJets"),
		minNumber = cms.uint32(2)
		)

ptEtaRestrictedSeq = cms.Sequence(
		ptEtaRestrictedRecoLeptons
		*ptEtaRestrictedRecoLeptonsFilter
		*ptEtaRestrictedLeadRecoLepton
		*ptEtaRestrictedLeadRecoLeptonFilter
		*ptEtaRestrictedRecoJets
		*ptEtaRestrictedRecoJetsFilter
		)

## end modules which apply pT and eta cuts to lepton and jets

## these modules apply the dilepton mass cut to matched reco leptons
recoDileptonCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("ptEtaRestrictedLeadRecoLepton ptEtaRestrictedRecoLeptons"),
		role = cms.string("leadingLepton subleadingLepton"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 200 && daughter(0).pt > daughter(1).pt")
		)
recoDileptonCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("recoDileptonCandidate"),
		minNumber = cms.uint32(1)
		)

recoDileptonCandidateSeq = cms.Sequence(recoDileptonCandidate*recoDileptonCandidateFilter)

## end modules which apply dilepton mass cut to reco leptons



