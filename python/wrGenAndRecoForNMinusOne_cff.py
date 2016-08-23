##based on ttBarGenAndRecoForNMinusOne_cff
import FWCore.ParameterSet.Config as cms
from ttBarGenAndRecoForNMinusOne_cff import *

##############################################################################################
##modules and sequences for selecting gen leptons and gen quarks with mother requirements

matchedGenEleFromWR = bareGenEle.clone(
		cut = cms.string("abs(pdgId) == 11 && ( (abs(mother(0).pdgId) == 23 && mother(0).status == 62 ) || status == 23 )")
		)
matchedGenEleFromWRFilter = bareGenEleFilter.clone(src = cms.InputTag("matchedGenEleFromWR"))

matchedGenPartonFromWR = bareGenEle.clone(cut = cms.string("(abs(pdgId) < 7 || abs(pdgId) == 21) && status == 23"))
matchedGenPartonFromWRFilter = bareGenEleFilter.clone(src = cms.InputTag("matchedGenPartonFromWR"))

matchedGenMuonFromDY = bareGenEle.clone(
		cut = cms.string("abs(pdgId) == 13 && ( (abs(mother(0).pdgId) == 23 && mother(0).status == 62 ) || status == 23 )")
		)
matchedGenMuonFromDYFilter = bareGenEleFilter.clone(src = cms.InputTag("matchedGenMuonFromDY"))

matchedGenEleAndQuarkFromDYSeq = cms.Sequence(
		matchedGenEleFromWR*matchedGenEleFromWRFilter
		*matchedGenPartonFromWR*matchedGenPartonFromWRFilter
		)

matchedGenMuAndQuarkFromDYSeq = cms.Sequence(
		matchedGenMuonFromDY*matchedGenMuonFromDYFilter
		*matchedGenPartonFromWR*matchedGenPartonFromWRFilter
		)
##end modules and sequences for selecting gen leptons and gen quarks with mother requirements



##############################################################################################
##modules and sequences for selecting gen leptons and gen jets with mother requirements
##the gen jets are dR matched to the gen quarks which come from the DY hard interaction
matchedGenJetFromDY = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedGenJetsDYNoCuts"),
		dRforMatching = cms.double(0.6),
		treeName = cms.string("matchedGenJetsDYNoCutsTree"),
		lowLevelCollTag = cms.InputTag("matchedGenPartonFromWR"),
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



