import FWCore.ParameterSet.Config as cms
from genElectronChannelModules_cff import *

printDYJetsParticleTree = printParticleTree.clone(
		src = cms.InputTag("prunedGenParticles")
		)

### begin modules which select GEN electrons, GEN jets, quarks, and gluons coming from the hard interaction or Z decay

#bareMatchedLeadingGenEle is an EDFilter module
#if the filter returns False, the path running the filter will stop
dyJetsBareMatchedGenEle = bareMatchedLeadingGenEle.clone(
		cut = cms.string("abs(pdgId) == 11 && abs(mother(0).pdgId) == 23 && mother(0).status == 62")
		)

dyJetsBareMatchedGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("dyJetsBareMatchedGenEle"),
		minNumber = cms.uint32(2)
		)

dyJetsBareGenJets = bareGenJet.clone()

dyJetsBareGenJetsFilter = dyJetsBareMatchedGenEleFilter.clone(
		src = cms.InputTag("dyJetsBareGenJets")
		)

dyJetsBareMatchedGenQuark = dyJetsBareMatchedGenEle.clone(
		cut = cms.string("abs(pdgId) < 7 && status == 23")
		)
#this filter is not used in the sequence below, but is used to study the gen quarks which leave the hard interaction
dyJetsBareMatchedGenQuarkFilter = dyJetsBareMatchedGenEleFilter.clone(
		src = cms.InputTag("dyJetsBareMatchedGenQuark")
		)

dyJetsBareMatchedGenGluon = dyJetsBareMatchedGenQuark.clone(
		cut = cms.string("abs(pdgId) == 21 && status == 23")
		)

#this filter is not used in the sequence below, but is used to study the gen gluons which leave the hard interaction
dyJetsBareMatchedGenGluonFilter = dyJetsBareMatchedGenEleFilter.clone(
		src = cms.InputTag("dyJetsBareMatchedGenGluon"),
		minNumber = cms.uint32(1)
		)

dyJetsBareMatchedGenLeptonsAndPartonsSeq = cms.Sequence(
		dyJetsBareMatchedGenEle*dyJetsBareMatchedGenEleFilter
		+dyJetsBareGenJets*dyJetsBareGenJetsFilter
		+dyJetsBareMatchedGenQuark
		+dyJetsBareMatchedGenGluon
		)

### end modules which select GEN electrons, GEN jets, quarks, and gluons coming from the hard interaction or Z decay

### begin modules which merge GEN quark and gluon collections into a single collection, and
### match objects in the gen jets collection to these gen quarks and gluons
#CandViewMerger will work if one of the input collections is empty
dyJetsMergeGenMatchedPartons = cms.EDProducer("CandViewMerger",
		src = cms.VInputTag("dyJetsBareMatchedGenQuark","dyJetsBareMatchedGenGluon")
		)

dyJetsCheckMergedGenMatchedPartons = cms.EDFilter("checkParticleContent",
		inputParticlesCollTag = cms.InputTag("dyJetsMergeGenMatchedPartons")
		)

dyJetsMatchGenJetsToGenPartons = matchGenJetsToGenQuarksNoCuts.clone(
		matchedOutputCollectionName = cms.string("matchedGenJets"),
		lowLevelCollTag = cms.InputTag("dyJetsMergeGenMatchedPartons"),
		higherLevelCollTag = cms.InputTag("dyJetsBareGenJets")
		)

dyJetsMatchGenJetsToGenPartonsFilter = matchGenJetsToGenQuarksNoCutsFilter.clone(
		src = cms.InputTag("dyJetsMatchGenJetsToGenPartons","matchedGenJets")
		)

dyJetsMatchGenJetsToGenPartonsSeq = cms.Sequence(
		dyJetsMergeGenMatchedPartons*dyJetsCheckMergedGenMatchedPartons
		*dyJetsMatchGenJetsToGenPartons*dyJetsMatchGenJetsToGenPartonsFilter
		)

### end modules which merge gen quarks with gen gluons, and match gen jets to the gen partons via dR

### begin modules which apply pT and eta cuts to the gen electrons and gen jets
dyJetsPtEtaRestrictedMatchedGenLeptons = simultaneousPtEtaCutMatchedSubleadingGenEle.clone(
		src = cms.InputTag("dyJetsBareMatchedGenEle")
		)

dyJetsPtEtaRestrictedMatchedGenLeptonsFilter = dyJetsBareMatchedGenEleFilter.clone(
		src = cms.InputTag("dyJetsPtEtaRestrictedMatchedGenLeptons")
		)

dyJetsPtEtaRestrictedMatchedGenJets = simultaneousPtEtaCutMatchedGenJets.clone(
		src = cms.InputTag("dyJetsMatchGenJetsToGenPartons","matchedGenJets")
		)

dyJetsPtEtaRestrictedMatchedGenJetsFilter = dyJetsBareGenJetsFilter.clone(
		src = cms.InputTag("dyJetsPtEtaRestrictedMatchedGenJets")
		)

dyJetsPtEtaRestrictedMatchedGenElesAndJetsSeq = cms.Sequence(
		dyJetsPtEtaRestrictedMatchedGenLeptons*dyJetsPtEtaRestrictedMatchedGenLeptonsFilter
		*dyJetsPtEtaRestrictedMatchedGenJets*dyJetsPtEtaRestrictedMatchedGenJetsFilter
		)

### end modules which apply pT and eta cuts to the gen electrons and gen jets

### begin modules which apply dR(lepton, jet) separation cut
dyJetsMatchedGenJetLeptonDrSeparation = cms.EDProducer("applyLeptonJetDrCut",
		outputLeptonsCollectionName = cms.string("leptonsPassingDrCut"),
		outputJetsCollectionName = cms.string("jetsPassingDrCut"),
		minDrSeparation = cms.double(0.4),
		minDileptonMassCut = cms.double(-1),
		inputJetsCollTag = cms.InputTag("dyJetsPtEtaRestrictedMatchedGenJets"),
		inputLeptonsCollTag = cms.InputTag("dyJetsPtEtaRestrictedMatchedGenLeptons")
		)

dyJetsMatchedGenJetLeptonDrSeparationGenJetsFilter = dyJetsBareMatchedGenEleFilter.clone(
		src = cms.InputTag("dyJetsMatchedGenJetLeptonDrSeparation","jetsPassingDrCut")
		)

dyJetsMatchedGenJetLeptonDrSeparationGenLeptonsFilter = dyJetsMatchedGenJetLeptonDrSeparationGenJetsFilter.clone(
		src = cms.InputTag("dyJetsMatchedGenJetLeptonDrSeparation","leptonsPassingDrCut")
		)

dyJetsMatchedGenJetLeptonDrSeparationSeq = cms.Sequence(
		dyJetsMatchedGenJetLeptonDrSeparation
		*dyJetsMatchedGenJetLeptonDrSeparationGenJetsFilter
		*dyJetsMatchedGenJetLeptonDrSeparationGenLeptonsFilter
		)
### end modules which apply dR(lepton, jet) separation cut


### begin modules which require one lepton with high pT
dyJetsRequireHighPtMatchedGenLepton = requireGenMatchedHighPtEle.clone(
		src = cms.InputTag("dyJetsMatchedGenJetLeptonDrSeparation","leptonsPassingDrCut")
		)

dyJetsRequireHighPtMatchedGenLeptonFilter = requireGenMatchedHighPtEleFilter.clone(src = cms.InputTag("dyJetsRequireHighPtMatchedGenLepton"))

dyJetsRequireHighPtMatchedGenLeptonSeq = cms.Sequence(
		dyJetsRequireHighPtMatchedGenLepton*dyJetsRequireHighPtMatchedGenLeptonFilter
		)

### end modules which require one lepton with high pT


### begin modules which select all leptons which pass dR(L,J) cut
dyJetsSelectLeptonsAfterDrSeparation = dyJetsBareGenJets.clone(src = cms.InputTag("dyJetsMatchedGenJetLeptonDrSeparation","leptonsPassingDrCut"))

dyJetsSelectLeptonsAfterDrSeparationFilter = dyJetsBareGenJetsFilter.clone(src = cms.InputTag("dyJetsSelectLeptonsAfterDrSeparation"))

dyJetsSelectLeptonsAfterDrSeparationSeq = cms.Sequence(
		dyJetsSelectLeptonsAfterDrSeparation*dyJetsSelectLeptonsAfterDrSeparationFilter
		)
### end modules which select all leptons which pass dR(L,J) cut



### begin modules which apply dilepton mass cut
dyJetsMatchedGenDileptonCand = genMatchedDiLeptonCandidate.clone(
		decay = cms.string("dyJetsRequireHighPtMatchedGenLepton dyJetsSelectLeptonsAfterDrSeparation")
		)

dyJetsMatchedGenDileptonCandFilter = genMatchedDiLeptonCandidateFilter.clone(src = cms.InputTag("dyJetsMatchedGenDileptonCand"))

dyJetsMatchedGenDileptonCandSeq = cms.Sequence(
		dyJetsMatchedGenDileptonCand*dyJetsMatchedGenDileptonCandFilter
		)

### end modules which apply dilepton mass cut


### begin modules which apply four object mass cut
dyJetsMatchedGenFourObjMass = cms.EDProducer("applyFourObjMassCut",
		outputJetsCollectionName = cms.string("jetsPassingFourObjMassCut"),
		outputLeptonsCollectionName = cms.string("leptonsPassingFourObjMassCut"),
		minFourObjMassCut = cms.double(600),
		minDileptonMassCut = cms.double(200),
		inputJetsCollTag = cms.InputTag("dyJetsMatchedGenJetLeptonDrSeparation","jetsPassingDrCut"),
		inputLeptonsCollTag = cms.InputTag("dyJetsMatchedGenJetLeptonDrSeparation","leptonsPassingDrCut")
		)

dyJetsMatchedGenFourObjMassLeptonsFilter = dyJetsBareMatchedGenEleFilter.clone(src = cms.InputTag("dyJetsMatchedGenFourObjMass","leptonsPassingFourObjMassCut"))

dyJetsMatchedGenFourObjMassJetsFilter = dyJetsMatchedGenFourObjMassLeptonsFilter.clone(src = cms.InputTag("dyJetsMatchedGenFourObjMass","jetsPassingFourObjMassCut"))

dyJetsMatchedGenFourObjMassSeq = cms.Sequence(
		dyJetsMatchedGenFourObjMass*dyJetsMatchedGenFourObjMassLeptonsFilter*dyJetsMatchedGenFourObjMassJetsFilter
		)

### end modules which apply four object mass cut




