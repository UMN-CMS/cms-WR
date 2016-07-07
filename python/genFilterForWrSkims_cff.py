import FWCore.ParameterSet.Config as cms

#################################################

#discard any events in which a Nu_Rtau appears in the decay chain   pdgId 9900016
hasGenNuTau = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) == 9900016")
		)

hasGenNuTauFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("hasGenNuTau"),
		minNumber = cms.uint32(1)
		)

skipGenNuTauSeq = cms.Sequence(hasGenNuTau *~hasGenNuTauFilter)

#discard any events in which a Nu_Rmu appears in the decay chain   pdgId 9900014
hasGenNuMu = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) == 9900014")
		)

hasGenNuMuFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("hasGenNuMu"),
		minNumber = cms.uint32(1)
		)

skipGenNuMuSeq = cms.Sequence(hasGenNuMu *~hasGenNuMuFilter)

#discard any events in which a Nu_Re appears in the decay chain   pdgId 9900012
hasGenNuEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) == 9900012")
		)

hasGenNuEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("hasGenNuEle"),
		minNumber = cms.uint32(1)
		)

skipGenNuEleSeq = cms.Sequence(hasGenNuEle *~hasGenNuEleFilter)
#################################################

#################################################
#WR producers and filters

#select the WR with status 62 that decays to the first lepton and heavy Nu in the event. In the centrally produced samples, this WR
#has one mother which is a WR produced by the hard interaction with status 22. The privately generated samples will not pass this filter.
trueWR = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) == 9900024 && abs(status) == 62 && abs(mother(0).pdgId) == 9900024 && abs(mother(0).status) == 22")
		)

trueWRFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("trueWR"),
		minNumber = cms.uint32(1)
		)

trueWRSeq = cms.Sequence(trueWR * trueWRFilter)
#################################################

#################################################
#Nu producers and filters

#select the Nu with status 22 which has WR mother with status 62, or the Nu with status 52 with Nu mother with status 22 in
#events where the lepton from the WR emits a photon
trueNuEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("(abs(pdgId) == 9900012 && abs(status) == 22 && abs(mother(0).pdgId) == 9900024) || (abs(pdgId) == 9900012 && abs(status) == 52 && abs(mother(0).pdgId) == 9900012 && abs(mother(0).status) == 22)")
		)

trueNuMu = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("(abs(pdgId) == 9900014 && abs(status) == 22 && abs(mother(0).pdgId) == 9900024) || (abs(pdgId) == 9900014 && abs(status) == 52 && abs(mother(0).pdgId) == 9900014 && abs(mother(0).status) == 22)")
		)

mergeTrueNuMuEle = cms.EDProducer("CandViewMerger",
		src = cms.VInputTag("trueNuEle","trueNuMu")
		)

mergeTrueNuMuEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("mergeTrueNuMuEle"),
		minNumber = cms.uint32(1)
		)


#skip the event if a nu_R is found with status 22 which is the daughter of a nu_R with status 22 or 52 and different abs(pdgId)
flavorChangeVeto = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) > 9900011 && abs(pdgId) < 9900020 && abs(status) == 22 && abs(mother(0).pdgId) > 9900011 && abs(mother(0).pdgId) < 9900020 && (abs(mother(0).status) == 22 || abs(mother(0).status) == 52) && abs(abs(pdgId) - abs(mother(0).pdgId)) > 0")
		)

flavorChangeVetoFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("flavorChangeVeto"),
		minNumber = cms.uint32(1)
		)

trueNuAndFlavorChangeVetoSeq = cms.Sequence(
		#trueNuEle + trueNuMu
		#+ mergeTrueNuMuEle * mergeTrueNuMuEleFilter
		flavorChangeVeto *~flavorChangeVetoFilter
		)
#################################################


#################################################
#filters and producers for quarks and lepton from Nu decay via virtual WR

#select events where the the lepton from the Nu decay has the same flavor as the Nu
trueLeptonFromNu = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("(abs(pdgId) == 11 || abs(pdgId) == 13) && (abs(mother(0).pdgId) - 9900001 - abs(pdgId)) == 0")
		)

trueLeptonFromNuFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("trueLeptonFromNu"),
		minNumber = cms.uint32(1)
		)

trueLeptonFromNuSeq = cms.Sequence(trueLeptonFromNu * trueLeptonFromNuFilter)

#select events where two quarks are produced from a Nu (requiring status 22 or 52 will not help)
trueQuarksFromNu = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) < 7 && abs(mother(0).pdgId) > 9900011")
		)

trueQuarksFromNuFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("trueQuarksFromNu"),
		minNumber = cms.uint32(2)
		)

trueQuarksFromNuSeq = cms.Sequence(trueQuarksFromNu * trueQuarksFromNuFilter)
#################################################


#################################################
#master sequence combining all sequences in this file

completeGenFilterSeq = cms.Sequence(
		skipGenNuTauSeq
		#*trueWRSeq
		*trueNuAndFlavorChangeVetoSeq
		*trueLeptonFromNuSeq
		*trueQuarksFromNuSeq
		)

#################################################




