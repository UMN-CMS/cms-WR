import FWCore.ParameterSet.Config as cms

##no cuts, just transforming the different collections into reco::Candidate collections
genParticlesFormatter = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("prunedGenParticles"),
		src = cms.InputTag("genParticles"),
		cut = cms.string("")
		)
genParticlesFormatterFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genParticlesFormatter"),
		minNumber = cms.uint32(2)
		)

recoLeptonsFormatter = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedElectrons"),
		#src = cms.InputTag("slimmedMuons"),
		cut = cms.string("")
		)
recoLeptonsFormatterFilter = cms.EDFilter("",
		src = cms.InputTag("recoLeptonsFormatter"),
		minNumber = cms.uint32(2)
		)

recoJetsFormatter = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedJets"),
		cut = cms.string("")
		)
recoJetsFormatterFilter = cms.EDFilter("",
		src = cms.InputTag("recoJetsFormatter"),
		minNumber = cms.uint32(2)
		)

genAndRecoFormatterSeq = cms.Sequence(
		genParticlesFormatter*genParticlesFormatterFilter
		#*recoLeptonsFormatter*recoLeptonsFormatterFilter
		#*recoJetsFormatter*recoJetsFormatterFilter
		)


##with pt and eta cuts on gen leptons and quarks
genParticlesPtEtaConstrainedLeptons = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("genParticlesFormatter"),
		cut = cms.string("abs(pdgId) == 11 && pt>40 && abs(eta) < 2.4")
		)
genParticlesPtEtaConstrainedLeptonsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genParticlesPtEtaConstrainedLeptons"),
		minNumber = cms.uint32(2)
		)

genParticlesPtEtaConstrainedQuarks = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("genParticlesFormatter"),
		cut = cms.string("abs(pdgId) < 7 && abs(pdgId) > 0 && pt>40 && abs(eta) < 2.5")
		)
genParticlesPtEtaConstrainedQuarksFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genParticlesPtEtaConstrainedQuarks"),
		minNumber = cms.uint32(2)
		)

genParticlesPtEtaConstrainedSeq = cms.Sequence(
		genParticlesPtEtaConstrainedLeptons*genParticlesPtEtaConstrainedLeptonsFilter
		*genParticlesPtEtaConstrainedQuarks*genParticlesPtEtaConstrainedQuarksFilter
		)



