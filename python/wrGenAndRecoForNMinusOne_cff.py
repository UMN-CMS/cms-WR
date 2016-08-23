##based on ttBarGenAndRecoForNMinusOne_cff
import FWCore.ParameterSet.Config as cms
from ttBarGenAndRecoForNMinusOne_cff import *

##############################################################################################
##modules and sequences for selecting gen leptons and gen quarks with mother requirements

matchedGenEleFromWR = bareGenEle.clone(
		cut = cms.string("abs(pdgId) == 11 && (abs(mother(0).pdgId) == 9900024 || abs(mother(0).pdgId) == 9900012) ")
		)
matchedGenEleFromWRFilter = bareGenEleFilter.clone(src = cms.InputTag("matchedGenEleFromWR"))

matchedGenPartonFromWREle = bareGenEle.clone(cut = cms.string("abs(pdgId) < 7 && abs(mother(0).pdgId) == 9900012"))
matchedGenPartonFromWREleFilter = bareGenEleFilter.clone(src = cms.InputTag("matchedGenPartonFromWREle"))

matchedGenMuonFromWR = bareGenEle.clone(
		cut = cms.string("abs(pdgId) == 13 && (abs(mother(0).pdgId) == 9900024 || abs(mother(0).pdgId) == 9900014)")
		)
matchedGenMuonFromWRFilter = bareGenEleFilter.clone(src = cms.InputTag("matchedGenMuonFromWR"))

matchedGenPartonFromWRMu = bareGenEle.clone(cut = cms.string("abs(pdgId) < 7 && abs(mother(0).pdgId) == 9900014"))
matchedGenPartonFromWRMuFilter = bareGenEleFilter.clone(src = cms.InputTag("matchedGenPartonFromWRMu"))

matchedGenEleAndQuarkFromWRSeq = cms.Sequence(
		matchedGenEleFromWR*matchedGenEleFromWRFilter
		*matchedGenPartonFromWREle*matchedGenPartonFromWREleFilter
		)

matchedGenMuAndQuarkFromWRSeq = cms.Sequence(
		matchedGenMuonFromWR*matchedGenMuonFromWRFilter
		*matchedGenPartonFromWRMu*matchedGenPartonFromWRMuFilter
		)
##end modules and sequences for selecting gen leptons and gen quarks with mother requirements



##############################################################################################
##modules and sequences for selecting gen leptons and gen jets with mother requirements
##the gen jets are dR matched to the gen quarks which come from the hard interaction
matchedGenJetFromWREle = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedGenJetsWREleNoCuts"),
		dRforMatching = cms.double(0.6),
		treeName = cms.string("matchedGenJetsWREleNoCutsTree"),
		lowLevelCollTag = cms.InputTag("matchedGenPartonFromWREle"),
		higherLevelCollTag = cms.InputTag("bareGenJet")
		)

matchedGenJetFromWREleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchedGenJetFromWREle","matchedGenJetsWREleNoCuts"),
		minNumber = cms.uint32(2)
		)

matchedGenJetWREleSeq = cms.Sequence(matchedGenJetFromWREle*matchedGenJetFromWREleFilter)

condenseMatchedGenJetWREleCollName = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("matchedGenJetFromWREle","matchedGenJetsWREleNoCuts"),
		cut = cms.string("")
		)
condenseMatchedGenJetWREleCollNameFilter = bareGenEleFilter.clone(src = cms.InputTag("condenseMatchedGenJetWREleCollName"))

condenseGenJetWREleSeq = cms.Sequence(condenseMatchedGenJetWREleCollName*condenseMatchedGenJetWREleCollNameFilter)


matchedGenJetFromWRMu = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedGenJetsWRMuNoCuts"),
		dRforMatching = cms.double(0.6),
		treeName = cms.string("matchedGenJetsWRMuNoCutsTree"),
		lowLevelCollTag = cms.InputTag("matchedGenPartonFromWRMu"),
		higherLevelCollTag = cms.InputTag("bareGenJet")
		)

matchedGenJetFromWRMuFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchedGenJetFromWRMu","matchedGenJetsWRMuNoCuts"),
		minNumber = cms.uint32(2)
		)

matchedGenJetWRMuSeq = cms.Sequence(matchedGenJetFromWRMu*matchedGenJetFromWRMuFilter)

condenseMatchedGenJetWRMuCollName = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("matchedGenJetFromWRMu","matchedGenJetsWRMuNoCuts"),
		cut = cms.string("")
		)
condenseMatchedGenJetWRMuCollNameFilter = bareGenEleFilter.clone(src = cms.InputTag("condenseMatchedGenJetWRMuCollName"))

condenseGenJetWRMuSeq = cms.Sequence(condenseMatchedGenJetWRMuCollName*condenseMatchedGenJetWRMuCollNameFilter)
##end modules and sequences for selecting gen leptons and gen jets with mother requirements



