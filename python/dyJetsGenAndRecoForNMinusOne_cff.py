##based on ttBarGenAndRecoForNMinusOne_cff
import FWCore.ParameterSet.Config as cms
from ttBarGenAndRecoForNMinusOne_cff import *

##############################################################################################
##modules and sequences for selecting gen leptons and gen jets without mother requirements
#bareGenEle = cms.EDFilter("CandViewSelector",
#		src = cms.InputTag("prunedGenParticles"),
#		cut = cms.string("abs(pdgId) == 11")
#		)
#
#bareGenEleFilter = cms.EDFilter("CandViewCountFilter",
#		src = cms.InputTag("bareGenEle"),
#		minNumber = cms.uint32(2)
#		)
#
#bareGenJet = cms.EDFilter("CandViewSelector",
#		src = cms.InputTag("slimmedGenJets"),
#		cut = cms.string("")
#		)
#
#bareGenJetFilter = cms.EDFilter("CandViewCountFilter",
#		src = cms.InputTag("bareGenJet"),
#		minNumber = cms.uint32(2)
#		)
#
#bareGenMu = bareGenEle.clone(
#		cut = cms.string("abs(pdgId) == 13")
#		)
#
#bareGenMuFilter = bareGenEleFilter.clone(
#		src = cms.InputTag("bareGenMu")
#		)
#
#genEleChnlGenJetsNoMatchingSequence = cms.Sequence(
#		bareGenEle*bareGenEleFilter
#		*bareGenJet*bareGenJetFilter
#		)
#
#genMuChnlGenJetsNoMatchingSequence = cms.Sequence(
#		bareGenMu*bareGenMuFilter
#		*bareGenJet*bareGenJetFilter
#		)
##end modules and sequence for selecting gen leptons and gen jets without mother requirements

###############################################################################################
###modules and sequences for selecting gen leptons and gen quarks without mother requirements
#bareGenQuark = bareGenEle.clone(cut = cms.string("abs(pdgId) < 7"))
#bareGenQuarkFilter = bareGenEleFilter.clone(src = cms.InputTag("bareGenQuark"))
#
#genEleChnlGenQuarksNoMatchingSequence = cms.Sequence(
#		bareGenEle*bareGenEleFilter
#		*bareGenQuark*bareGenQuarkFilter
#		)
#
#genMuChnlGenQuarksNoMatchingSequence = cms.Sequence(
#		bareGenMu*bareGenMuFilter
#		*bareGenQuark*bareGenQuarkFilter
#		)
###end modules and sequences for selecting gen leptons and gen quarks without mother requirements


##############################################################################################
##modules and sequences for selecting gen leptons and gen quarks with mother requirements

matchedGenEleFromDY = bareGenEle.clone(
		cut = cms.string("abs(pdgId) == 11 && ( (abs(mother(0).pdgId) == 23 && mother(0).status == 62 ) || status == 23 )")
		)
matchedGenEleFromDYFilter = bareGenEleFilter.clone(src = cms.InputTag("matchedGenEleFromDY"))

matchedGenPartonFromDY = bareGenEle.clone(cut = cms.string("(abs(pdgId) < 7 || abs(pdgId) == 21) && status == 23"))
#matchedGenPartonFromDYFilter = bareGenEleFilter.clone(src = cms.InputTag("matchedGenPartonFromDY"))
#######temporary test
matchedGenPartonFromDYFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchedGenPartonFromDY"),
		minNumber = cms.uint32(0)
		)
######

matchedGenMuonFromDY = bareGenEle.clone(
		cut = cms.string("abs(pdgId) == 13 && ( (abs(mother(0).pdgId) == 23 && mother(0).status == 62 ) || status == 23 )")
		)
matchedGenMuonFromDYFilter = bareGenEleFilter.clone(src = cms.InputTag("matchedGenMuonFromDY"))

matchedGenEleAndQuarkFromDYSeq = cms.Sequence(
		matchedGenEleFromDY*matchedGenEleFromDYFilter
		*matchedGenPartonFromDY*matchedGenPartonFromDYFilter
		)

matchedGenMuAndQuarkFromDYSeq = cms.Sequence(
		matchedGenMuonFromDY*matchedGenMuonFromDYFilter
		*matchedGenPartonFromDY*matchedGenPartonFromDYFilter
		)
##end modules and sequences for selecting gen leptons and gen quarks with mother requirements



##############################################################################################
##modules and sequences for selecting gen leptons and gen jets with mother requirements
##the gen jets are dR matched to the gen quarks which come from the DY hard interaction
matchedGenJetFromDY = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedGenJetsDYNoCuts"),
		dRforMatching = cms.double(0.6),
		treeName = cms.string("matchedGenJetsDYNoCutsTree"),
		lowLevelCollTag = cms.InputTag("matchedGenPartonFromDY"),
		higherLevelCollTag = cms.InputTag("bareGenJet")
		)

#matchedGenJetFromDYFilter = cms.EDFilter("CandViewCountFilter",
#		src = cms.InputTag("matchedGenJetFromDY","matchedGenJetsDYNoCuts"),
#		minNumber = cms.uint32(2)
#		)

#######temporary test
matchedGenJetFromDYFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchedGenJetFromDY","matchedGenJetsDYNoCuts"),
		minNumber = cms.uint32(0)
		)
#######

matchedGenJetDYSeq = cms.Sequence(matchedGenJetFromDY*matchedGenJetFromDYFilter)

condenseMatchedGenJetDYCollName = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("matchedGenJetFromDY","matchedGenJetsDYNoCuts"),
		cut = cms.string("")
		)
#condenseMatchedGenJetDYCollNameFilter = bareGenEleFilter.clone(src = cms.InputTag("condenseMatchedGenJetDYCollName"))

#######temporary test
condenseMatchedGenJetDYCollNameFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("condenseMatchedGenJetDYCollName"),
		minNumber = cms.uint32(0)
		)
#######

condenseGenJetDYSeq = cms.Sequence(condenseMatchedGenJetDYCollName*condenseMatchedGenJetDYCollNameFilter)
##end modules and sequences for selecting gen leptons and gen jets with mother requirements



###############################################################################################
###modules and sequences for selecting reco leptons and reco jets passing ID
#bareRecoJetPassingId = cms.EDFilter("CandViewSelector",
#		src = cms.InputTag("slimmedJets"),
#		cut = cms.string("(neutralHadronEnergyFraction<0.90 && neutralEmEnergyFraction<0.9 && (chargedMultiplicity+neutralMultiplicity)>1 && muonEnergyFraction<0.8) && ((abs(eta)<=2.4 && chargedHadronEnergyFraction>0 && chargedMultiplicity>0 && chargedEmEnergyFraction<0.90) || abs(eta)>2.4) ")
#		)
#
#bareRecoJetPassingIdFilter = bareGenEleFilter.clone(
#		src = cms.InputTag("bareRecoJetPassingId")
#		)
#
#heepRecoEle = cms.EDFilter("CandViewSelector",
#		src = cms.InputTag("HEEPIDSelector"),
#		cut = cms.string("")
#		)
#heepRecoEleFilter = bareGenEleFilter.clone(src = cms.InputTag("heepRecoEle"))
#
#recoEleJetPassingIdSeq = cms.Sequence(
#		bareRecoJetPassingId*bareRecoJetPassingIdFilter
#		*heepRecoEle*heepRecoEleFilter
#		)
#
#wrTunePMuProd = cms.EDProducer("TunePMuonProducer",
#		src = cms.InputTag("slimmedMuons")
#		)
#
#wrTunePMuProdFilter = cms.EDFilter("CandViewCountFilter",
#		src = cms.InputTag("wrTunePMuProd"),
#		minNumber = cms.uint32(2)
#		)
#
#wrTunePMuProdSeq = cms.Sequence(
#		wrTunePMuProd
#		*wrTunePMuProdFilter
#		)
#
## make a collection of TuneP muons which pass isHighPt ID
#isHighPtMuProd = cms.EDProducer("produceIsHighPtMuons",
#		src = cms.InputTag("wrTunePMuProd"),
#		outputCollectionName = cms.string("TunePMuonsPassingIsHighPtID"),
#		minPt = cms.double(0.)
#		)
#
#isHighPtMuProdFilter = cms.EDFilter("CandViewCountFilter",
#		src = cms.InputTag("isHighPtMuProd","TunePMuonsPassingIsHighPtID"),
#		minNumber = cms.uint32(2)
#		)
#
#isHighPtMuSeq = cms.Sequence(isHighPtMuProd*isHighPtMuProdFilter)
#
#highPtIdRecoMu = cms.EDFilter("CandViewSelector",
#		src = cms.InputTag("isHighPtMuProd","TunePMuonsPassingIsHighPtID"),
#		cut = cms.string("")
#		)
#highPtIdRecoMuFilter = heepRecoEleFilter.clone(src = cms.InputTag("highPtIdRecoMu"))
#
#recoMuJetPassingIdSeq = cms.Sequence(
#		bareRecoJetPassingId*bareRecoJetPassingIdFilter
#		*wrTunePMuProdSeq*isHighPtMuSeq
#		*highPtIdRecoMu*highPtIdRecoMuFilter
#		)
#
###end modules and sequences for selecting reco leptons and reco jets which pass ID


