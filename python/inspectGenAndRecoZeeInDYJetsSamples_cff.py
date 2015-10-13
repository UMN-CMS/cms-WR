import FWCore.ParameterSet.Config as cms

## select the two gen electrons with Z boson mother
genEleFromZ = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) == 11 && abs(mother(0).pdgId) == 23")
		)

genEleFromZFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genEleFromZ"),
		minNumber = cms.uint32(2)
		)

genZeeSeq = cms.Sequence(genEleFromZ * genEleFromZFilter)

## end modules to select the two gen electrons with Z boson mother 

## select the two reco electrons which are dR matched to the two gen electrons
## with Z boson mother
selectHEEPEles = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("HEEPIDSelector"),
		cut = cms.string("")
		)

selectHEEPElesFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("selectHEEPEles"),
		minNumber = cms.uint32(2)
		)

selectHEEPElesSeq = cms.Sequence(selectHEEPEles * selectHEEPElesFilter)

matchedRecoElesToGenZeeEles = cms.EDProducer("FindHigherLevelMatchedObject",
		matchedOutputCollectionName = cms.string("matchedRecoElesFromZee"),
		dRforMatching = cms.double(0.2),
		treeName = cms.string("matchingParametersRecoZedElesTree"),
		lowLevelCollTag = cms.InputTag("genEleFromZ"),
		higherLevelCollTag = cms.InputTag("selectHEEPEles")
		)

matchedRecoElesToGenZeeElesFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchedRecoElesToGenZeeEles","matchedRecoElesFromZee"),
		minNumber = cms.uint32(2)
		)

matchedRecoElesSeq = cms.Sequence(matchedRecoElesToGenZeeEles * matchedRecoElesToGenZeeElesFilter)
## end modules to select the two reco electrons which are dR matched to the two
## gen electrons with Z boson mother
