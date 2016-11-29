##based on ttBarGenAndRecoForNMinusOne_cff
import FWCore.ParameterSet.Config as cms
from ttBarGenAndRecoForNMinusOne_cff import *

##############################################################################################
##modules and sequences for selecting gen leptons and gen quarks with mother requirements

matchedGenEleFromDY = bareGenEle.clone(
		cut = cms.string("abs(pdgId) == 11 && ( (abs(mother(0).pdgId) == 23 && mother(0).status == 62 ) || status == 23 )")
		)
matchedGenEleFromDYFilter = bareGenEleFilter.clone(src = cms.InputTag("matchedGenEleFromDY"))

matchedGenPartonFromDY = bareGenEle.clone(cut = cms.string("(abs(pdgId) < 7 || abs(pdgId) == 21) && status == 23"))
matchedGenPartonFromDYFilter = bareGenEleFilter.clone(src = cms.InputTag("matchedGenPartonFromDY"))

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

matchedGenJetFromDYFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchedGenJetFromDY","matchedGenJetsDYNoCuts"),
		minNumber = cms.uint32(2)
		)


matchedGenJetDYSeq = cms.Sequence(matchedGenJetFromDY*matchedGenJetFromDYFilter)

condenseMatchedGenJetDYCollName = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("matchedGenJetFromDY","matchedGenJetsDYNoCuts"),
		cut = cms.string("")
		)
condenseMatchedGenJetDYCollNameFilter = bareGenEleFilter.clone(src = cms.InputTag("condenseMatchedGenJetDYCollName"))


condenseGenJetDYSeq = cms.Sequence(condenseMatchedGenJetDYCollName*condenseMatchedGenJetDYCollNameFilter)
##end modules and sequences for selecting gen leptons and gen jets with mother requirements



