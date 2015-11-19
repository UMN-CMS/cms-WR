import FWCore.ParameterSet.Config as cms
from genElectronChannelModules_cff import *

printTTBarParticleTree = printParticleTree.clone(
		src = cms.InputTag("prunedGenParticles")
		)

### begin modules which select GEN electrons, GEN jets, quarks, and gluons coming from the TTBar decay
ttBarBareMatchedGenW = bareMatchedLeadingGenEle.clone(
		cut = cms.string("abs(pdgId) == 24 && abs(mother(0).pdgId) == 6 && mother(0).status == 62")
		)

ttBarBareMatchedGenWFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ttBarBareMatchedGenW"),
		minNumber = cms.uint32(2)
		)

ttBarBareMatchedGenEle = bareMatchedLeadingGenEle.clone(
		cut = cms.string("abs(pdgId) == 11 && abs(mother(0).pdgId) == 24 && (mother(0).status == 52 || mother(0).status == 22)")
		)

ttBarBareMatchedGenEleFilter = ttBarBareMatchedGenWFilter.clone(
		src = cms.InputTag("ttBarBareMatchedGenEle")
		)

ttBarBareGenJets = bareGenJet.clone()

ttBarBareGenJetsFilter = ttBarBareMatchedGenEleFilter.clone(
		src = cms.InputTag("ttBarBareGenJets")
		)

ttBarBareMatchedGenQuark = ttBarBareMatchedGenEle.clone(
		cut = cms.string("abs(pdgId) == 5 && status == 23 && abs(mother(0).pdgId) == 6 && mother(0).status == 62")
		)

ttBarBareMatchedGenQuarkFilter = ttBarBareMatchedGenEleFilter.clone(
		src = cms.InputTag("ttBarBareMatchedGenQuark")
		)

ttBarBareMatchedGenLeptonsAndPartonsSeq = cms.Sequence(
		ttBarBareMatchedGenW*ttBarBareMatchedGenWFilter
		*ttBarBareMatchedGenEle*ttBarBareMatchedGenEleFilter
		+ttBarBareGenJets*ttBarBareGenJetsFilter
		+ttBarBareMatchedGenQuark*ttBarBareMatchedGenQuarkFilter
		)

### end modules which select GEN electrons, GEN jets, and quarks coming from the TTBar decay


### begin modules which match objects in the gen jets collection to gen quarks (b quarks)
ttBarMatchGenJetsToGenPartons = matchGenJetsToGenQuarksNoCuts.clone(
		matchedOutputCollectionName = cms.string("matchedGenJets"),
		lowLevelCollTag = cms.InputTag("ttBarBareMatchedGenQuark"),
		higherLevelCollTag = cms.InputTag("ttBarBareGenJets")
		)

ttBarMatchGenJetsToGenPartonsFilter = matchGenJetsToGenQuarksNoCutsFilter.clone(
		src = cms.InputTag("ttBarMatchGenJetsToGenPartons","matchedGenJets")
		)

ttBarMatchGenJetsToGenPartonsSeq = cms.Sequence(
		ttBarMatchGenJetsToGenPartons*ttBarMatchGenJetsToGenPartonsFilter
		)

### end modules which merge gen quarks with gen gluons, and match gen jets to the gen partons via dR


### begin modules which apply pT and eta cuts to the gen electrons and gen jets
ttBarPtEtaRestrictedMatchedGenLeptons = simultaneousPtEtaCutMatchedSubleadingGenEle.clone(
		src = cms.InputTag("ttBarBareMatchedGenEle")
		)

ttBarPtEtaRestrictedMatchedGenLeptonsFilter = ttBarBareMatchedGenEleFilter.clone(
		src = cms.InputTag("ttBarPtEtaRestrictedMatchedGenLeptons")
		)

ttBarPtEtaRestrictedMatchedGenJets = simultaneousPtEtaCutMatchedGenJets.clone(
		src = cms.InputTag("ttBarMatchGenJetsToGenPartons","matchedGenJets")
		)

ttBarPtEtaRestrictedMatchedGenJetsFilter = ttBarBareGenJetsFilter.clone(
		src = cms.InputTag("ttBarPtEtaRestrictedMatchedGenJets")
		)

ttBarPtEtaRestrictedMatchedGenElesAndJetsSeq = cms.Sequence(
		ttBarPtEtaRestrictedMatchedGenLeptons*ttBarPtEtaRestrictedMatchedGenLeptonsFilter
		*ttBarPtEtaRestrictedMatchedGenJets*ttBarPtEtaRestrictedMatchedGenJetsFilter
		)

### end modules which apply pT and eta cuts to the gen electrons and gen jets

### begin modules which apply dR(lepton, jet) separation cut
ttBarMatchedGenJetLeptonDrSeparation = cms.EDProducer("applyLeptonJetDrCut",
		outputLeptonsCollectionName = cms.string("leptonsPassingDrCut"),
		outputJetsCollectionName = cms.string("jetsPassingDrCut"),
		minDrSeparation = cms.double(0.4),
		minDileptonMassCut = cms.double(-1),
		inputJetsCollTag = cms.InputTag("ttBarPtEtaRestrictedMatchedGenJets"),
		inputLeptonsCollTag = cms.InputTag("ttBarPtEtaRestrictedMatchedGenLeptons")
		)

ttBarMatchedGenJetLeptonDrSeparationGenJetsFilter = ttBarBareMatchedGenEleFilter.clone(
		src = cms.InputTag("ttBarMatchedGenJetLeptonDrSeparation","jetsPassingDrCut")
		)

ttBarMatchedGenJetLeptonDrSeparationGenLeptonsFilter = ttBarMatchedGenJetLeptonDrSeparationGenJetsFilter.clone(
		src = cms.InputTag("ttBarMatchedGenJetLeptonDrSeparation","leptonsPassingDrCut")
		)

ttBarMatchedGenJetLeptonDrSeparationSeq = cms.Sequence(
		ttBarMatchedGenJetLeptonDrSeparation
		*ttBarMatchedGenJetLeptonDrSeparationGenJetsFilter
		*ttBarMatchedGenJetLeptonDrSeparationGenLeptonsFilter
		)
### end modules which apply dR(lepton, jet) separation cut


### begin modules which require one lepton with high pT
ttBarRequireHighPtMatchedGenLepton = requireGenMatchedHighPtEle.clone(
		src = cms.InputTag("ttBarMatchedGenJetLeptonDrSeparation","leptonsPassingDrCut")
		)

ttBarRequireHighPtMatchedGenLeptonFilter = requireGenMatchedHighPtEleFilter.clone(src = cms.InputTag("ttBarRequireHighPtMatchedGenLepton"))

ttBarRequireHighPtMatchedGenLeptonSeq = cms.Sequence(
		ttBarRequireHighPtMatchedGenLepton*ttBarRequireHighPtMatchedGenLeptonFilter
		)

### end modules which require one lepton with high pT


### begin modules which select all leptons which pass dR(L,J) cut
ttBarSelectLeptonsAfterDrSeparation = ttBarBareGenJets.clone(src = cms.InputTag("ttBarMatchedGenJetLeptonDrSeparation","leptonsPassingDrCut"))

ttBarSelectLeptonsAfterDrSeparationFilter = ttBarBareGenJetsFilter.clone(src = cms.InputTag("ttBarSelectLeptonsAfterDrSeparation"))

ttBarSelectLeptonsAfterDrSeparationSeq = cms.Sequence(
		ttBarSelectLeptonsAfterDrSeparation*ttBarSelectLeptonsAfterDrSeparationFilter
		)
### end modules which select all leptons which pass dR(L,J) cut



### begin modules which apply dilepton mass cut
ttBarMatchedGenDileptonCand = genMatchedDiLeptonCandidate.clone(
		decay = cms.string("ttBarRequireHighPtMatchedGenLepton ttBarSelectLeptonsAfterDrSeparation")
		)

ttBarMatchedGenDileptonCandFilter = genMatchedDiLeptonCandidateFilter.clone(src = cms.InputTag("ttBarMatchedGenDileptonCand"))

ttBarMatchedGenDileptonCandSeq = cms.Sequence(
		ttBarMatchedGenDileptonCand*ttBarMatchedGenDileptonCandFilter
		)

### end modules which apply dilepton mass cut


### begin modules which apply four object mass cut
ttBarMatchedGenFourObjMass = cms.EDProducer("applyFourObjMassCut",
		outputJetsCollectionName = cms.string("jetsPassingFourObjMassCut"),
		outputLeptonsCollectionName = cms.string("leptonsPassingFourObjMassCut"),
		minFourObjMassCut = cms.double(600),
		minDileptonMassCut = cms.double(200),
		inputJetsCollTag = cms.InputTag("ttBarMatchedGenJetLeptonDrSeparation","jetsPassingDrCut"),
		inputLeptonsCollTag = cms.InputTag("ttBarMatchedGenJetLeptonDrSeparation","leptonsPassingDrCut")
		)

ttBarMatchedGenFourObjMassLeptonsFilter = ttBarBareMatchedGenEleFilter.clone(src = cms.InputTag("ttBarMatchedGenFourObjMass","leptonsPassingFourObjMassCut"))

ttBarMatchedGenFourObjMassJetsFilter = ttBarMatchedGenFourObjMassLeptonsFilter.clone(src = cms.InputTag("ttBarMatchedGenFourObjMass","jetsPassingFourObjMassCut"))

ttBarMatchedGenFourObjMassSeq = cms.Sequence(
		ttBarMatchedGenFourObjMass*ttBarMatchedGenFourObjMassLeptonsFilter*ttBarMatchedGenFourObjMassJetsFilter
		)

### end modules which apply four object mass cut




