import FWCore.ParameterSet.Config as cms

##filters to select gen particles matched to WR decay products via pdgId

#this filter looks for muons whose mother is a WR
bareMatchedLeadingGenMuon = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) == 13 && abs(mother(0).pdgId) == 9900024")
		)

bareMatchedLeadingGenMuonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareMatchedLeadingGenMuon"),
		minNumber = cms.uint32(1)
		)

#this filter looks for muons whose mother is a heavy Nu
bareMatchedSubleadingGenMuon = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) == 13 && abs(mother(0).pdgId) == 9900014")
		)

bareMatchedSubleadingGenMuonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareMatchedSubleadingGenMuon"),
		minNumber = cms.uint32(1)
		)

#this filter looks for quarks whose real mother is a heavy Nu (virtuals are not tracked in Pythia) 
muBareMatchedGenQuark = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) < 7 && abs(mother(0).pdgId) == 9900014")
		)

muBareMatchedGenQuarkFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("muBareMatchedGenQuark"),
		minNumber = cms.uint32(2)
		)

muBareMatchedGenParticleSeq = cms.Sequence(bareMatchedLeadingGenMuon*bareMatchedLeadingGenMuonFilter*bareMatchedSubleadingGenMuon*
		bareMatchedSubleadingGenMuonFilter*muBareMatchedGenQuark*muBareMatchedGenQuarkFilter)


##filters to cut out gen particles which could never hit the tracker based on eta trajectory
##use these modules to determine efficiency of signal skim (skim applies pT cuts to leptons and possibly jets)
##skim efficiency on signal = # evts passing reco cuts / # evts before reco cuts with gen objects falling within tracker coverage
muEtaRestrictedMatchedGenQuark = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("muBareMatchedGenQuark"),
		cut = cms.string("abs(eta) < 2.5")
		)
muEtaRestrictedMatchedGenQuarkFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("muEtaRestrictedMatchedGenQuark"),
		minNumber = cms.uint32(2)
		)

etaRestrictedMatchedGenLeadingMuon = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareMatchedLeadingGenMuon"),
		cut = cms.string("abs(eta) < 2.4")
		)
etaRestrictedMatchedGenLeadingMuonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("etaRestrictedMatchedGenLeadingMuon"),
		minNumber = cms.uint32(1)
		)

etaRestrictedMatchedGenSubleadingMuon = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareMatchedSubleadingGenMuon"),
		cut = cms.string("abs(eta) < 2.4")
		)
etaRestrictedMatchedGenSubleadingMuonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("etaRestrictedMatchedGenSubleadingMuon"),
		minNumber = cms.uint32(1)
		)

muEtaRestrictedMatchedGenParticleSeq = cms.Sequence(
		muEtaRestrictedMatchedGenQuark
		*muEtaRestrictedMatchedGenQuarkFilter
		*etaRestrictedMatchedGenLeadingMuon
		*etaRestrictedMatchedGenLeadingMuonFilter
		*etaRestrictedMatchedGenSubleadingMuon
		*etaRestrictedMatchedGenSubleadingMuonFilter
		)


