import FWCore.ParameterSet.Config as cms

printParticleTree = cms.EDAnalyzer("ParticleListDrawer",
		maxEventsToPrint = cms.untracked.int32(-1),
		printVertex = cms.untracked.bool(False),
		printOnlyHardInteraction = cms.untracked.bool(False),
		#src = cms.InputTag("genParticles")
		src = cms.InputTag("prunedGenParticles")
		)

hasGenMuOrTau = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("(abs(pdgId) == 13 && abs(mother(0).pdgId) > 9900011) || (abs(pdgId) == 15 && abs(mother(0).pdgId) > 9900011)")
		)

hasGenMuOrTauFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("hasGenMuOrTau"),
		minNumber = cms.uint32(1)
		)

hasGenNuMuOrTau = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		#cut = cms.string("abs(pdgId) > 9900012 && abs(pdgId) < 9900024")
		cut = cms.string("abs(pdgId) == 9900012")
		)

hasGenNuMuOrTauFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("hasGenNuMuOrTau"),
		minNumber = cms.uint32(1)
		)

hasGenMuOrTauFlavorsSeq = cms.Sequence(
		#hasGenNuMuOrTau
		#*hasGenNuMuOrTauFilter
		hasGenMuOrTau*hasGenMuOrTauFilter
		)

bareGenJet = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedGenJets"),
		#src = cms.InputTag("ak4GenJets"),
		cut = cms.string("")
		)

bareGenJetFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareGenJet"),
		minNumber = cms.uint32(2)
		)

bareGenEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		#src = cms.InputTag("genParticles"),
		cut = cms.string("abs(pdgId) == 11")
		)

bareGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareGenEle"),
		minNumber = cms.uint32(2)
		)

bareGenParticleSeq = cms.Sequence(bareGenJet*bareGenJetFilter*bareGenEle*bareGenEleFilter)

bareGenJetSeq = cms.Sequence(bareGenJet*bareGenJetFilter)

## producers to find GenJets matched to Gen quarks
matchGenJetsToGenQuarksNoCuts = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedGenJetsNoCuts"),
		dRforMatching = cms.double(0.6),
		treeName = cms.string("matchedGenJetsNoCutsTree"),
		lowLevelCollTag = cms.InputTag("bareMatchedGenQuark"),
		higherLevelCollTag = cms.InputTag("bareGenJet")
		)

matchGenJetsToGenQuarksNoCutsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchGenJetsToGenQuarksNoCuts","matchedGenJetsNoCuts"),
		minNumber = cms.uint32(2)
		)

matchGenJetsToGenQuarksNoCutsNewPath = cms.EDProducer('FindHigherLevelMatchedObject',
		matchedOutputCollectionName = cms.string("matchedGenJetsNoCuts"),
		dRforMatching = cms.double(0.6),
		treeName = cms.string("matchedGenJetsNoCutsTree"),
		lowLevelCollTag = cms.InputTag("bareMatchedGenQuark"),
		higherLevelCollTag = cms.InputTag("bareGenJet")
		)

matchGenJetsToGenQuarksNoCutsNewPathFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("matchGenJetsToGenQuarksNoCutsNewPath","matchedGenJetsNoCuts"),
		minNumber = cms.uint32(2)
		)

matchGenJetsToGenQuarksSeq = cms.Sequence(matchGenJetsToGenQuarksNoCuts*matchGenJetsToGenQuarksNoCutsFilter)

## filters to select the generated Nu
bareMatchedNu = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("prunedGenParticles"),
		#src = cms.InputTag("genParticles"),
		cut = cms.string("abs(pdgId) == 9900012 && abs(mother(0).pdgId) == 9900024")
		)

bareMatchedNuFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareMatchedNu"),
		minNumber = cms.uint32(1)
		)

bareMatchedNuSeq = cms.Sequence(bareMatchedNu*bareMatchedNuFilter)


## filters to select the generated WR
bareMatchedWR = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("prunedGenParticles"),
		src = cms.InputTag("genParticles"),
		cut = cms.string("abs(pdgId) == 9900024 && abs(mother(0).pdgId) == 9900024 && abs(status) == 62")
		)

bareMatchedWRFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareMatchedWR"),
		minNumber = cms.uint32(1)
		)

bareMatchedWRSeq = cms.Sequence(bareMatchedWR*bareMatchedWRFilter)


##filters to select gen particles matched to WR decay products via pdgId
##heavyNu pdgId = 9900012, WR pdgId = 9900024

#this filter looks for electrons whose mother is a WR
bareMatchedLeadingGenEle = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("genParticles"),
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) == 11 && abs(mother(0).pdgId) == 9900024")
		)

bareMatchedLeadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareMatchedLeadingGenEle"),
		minNumber = cms.uint32(1)
		)

#this filter looks for electrons whose mother is a heavy Nu
bareMatchedSubleadingGenEle = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("genParticles"),
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) == 11 && abs(mother(0).pdgId) == 9900012")
		)

bareMatchedSubleadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareMatchedSubleadingGenEle"),
		minNumber = cms.uint32(1)
		)

#this filter looks for quarks whose real mother is a heavy Nu (virtuals are not tracked in Pythia) 
bareMatchedGenQuark = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("genParticles"),
		src = cms.InputTag("prunedGenParticles"),
		cut = cms.string("abs(pdgId) < 7 && abs(mother(0).pdgId) == 9900012")
		)

bareMatchedGenQuarkFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareMatchedGenQuark"),
		minNumber = cms.uint32(2)
		)

bareMatchedGenParticleSeq = cms.Sequence(
		bareMatchedLeadingGenEle
		*bareMatchedLeadingGenEleFilter
		*bareMatchedSubleadingGenEle
		*bareMatchedSubleadingGenEleFilter
		*bareMatchedGenQuark
		*bareMatchedGenQuarkFilter
		)


##filters to cut out gen particles which could never hit the tracker based on eta trajectory
##use these modules to determine efficiency of signal skim (skim applies pT cuts to leptons and possibly jets)
##skim efficiency on signal = # evts passing reco cuts / # evts before reco cuts with gen objects falling within tracker coverage

etaRestrictedMatchedGenQuark = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareMatchedGenQuark"),
		cut = cms.string("abs(eta) < 2.5")
		)
etaRestrictedMatchedGenQuarkFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("etaRestrictedMatchedGenQuark"),
		minNumber = cms.uint32(2)
		)

etaRestrictedMatchedGenLeadingLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareMatchedLeadingGenEle"),
		cut = cms.string("abs(eta) < 2.5")
		)
etaRestrictedMatchedGenLeadingLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("etaRestrictedMatchedGenLeadingLepton"),
		minNumber = cms.uint32(1)
		)

etaRestrictedMatchedGenSubleadingLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareMatchedSubleadingGenEle"),
		cut = cms.string("abs(eta) < 2.5")
		)
etaRestrictedMatchedGenSubleadingLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("etaRestrictedMatchedGenSubleadingLepton"),
		minNumber = cms.uint32(1)
		)

etaRestrictedMatchedGenParticleSeq = cms.Sequence(
		etaRestrictedMatchedGenQuark
		*etaRestrictedMatchedGenQuarkFilter
		*etaRestrictedMatchedGenLeadingLepton
		*etaRestrictedMatchedGenLeadingLeptonFilter
		*etaRestrictedMatchedGenSubleadingLepton
		*etaRestrictedMatchedGenSubleadingLeptonFilter
		)


## filters on pt and eta of gen leptons from WR decay, and gen jets
## matched to gen quarks from WR decay
simultaneousPtEtaCutMatchedLeadingGenEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareMatchedLeadingGenEle"),
		cut = cms.string("pt>40 && abs(eta) < 2.5")
		)

simultaneousPtEtaCutMatchedLeadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("simultaneousPtEtaCutMatchedLeadingGenEle"),
		minNumber = cms.uint32(1)
		)

simultaneousPtEtaCutMatchedSubleadingGenEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareMatchedSubleadingGenEle"),
		cut = cms.string("pt>40 && abs(eta) < 2.5")
		)

simultaneousPtEtaCutMatchedSubleadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("simultaneousPtEtaCutMatchedSubleadingGenEle"),
		minNumber = cms.uint32(1)
		)

simultaneousPtEtaCutMatchedGenJets = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("matchGenJetsToGenQuarksNoCuts","matchedGenJetsNoCuts"),
		cut = cms.string("pt>40 && abs(eta) < 2.5")
		)

simultaneousPtEtaCutMatchedGenJetsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("simultaneousPtEtaCutMatchedGenJets"),
		minNumber = cms.uint32(2)
		)


simultaneousPtEtaCutMatchedObjectsSeq = cms.Sequence(
		simultaneousPtEtaCutMatchedLeadingGenEle
		*simultaneousPtEtaCutMatchedLeadingGenEleFilter
		*simultaneousPtEtaCutMatchedSubleadingGenEle
		*simultaneousPtEtaCutMatchedSubleadingGenEleFilter
		*simultaneousPtEtaCutMatchedGenJets
		*simultaneousPtEtaCutMatchedGenJetsFilter
		)

## end simultaneous pt and eta filters on gen leptons from WR decays,
## and gen jets matched to gen quarks from WR decay


## apply the dR(lepton, jet) cut on gen leptons and gen jets matched to WR decay products
genMatchedJetLeptonDrSeparation = cms.EDProducer("applyLeptonJetDrCutMixedLeptonFlavor",
		outputJetsCollectionName = cms.string("genJetsPassingDrSeparationCut"),
		outputLeptonsOneCollectionName = cms.string("matchedLeadingGenElePassingDrSeparationCut"),
		outputLeptonsTwoCollectionName = cms.string("matchedSubleadingGenElePassingDrSeparationCut"),
		minDrSeparation = cms.double(0.4),
		inputJetsCollTag = cms.InputTag("simultaneousPtEtaCutMatchedGenJets"),
		inputLeptonsOneCollTag = cms.InputTag("simultaneousPtEtaCutMatchedLeadingGenEle"),
		inputLeptonsTwoCollTag = cms.InputTag("simultaneousPtEtaCutMatchedSubleadingGenEle"),
		minDileptonMassCut = cms.double(-1),
		minLeptonDrSeparation = cms.double(0.4)  #separation btwn the two leptons
		)

genMatchedJetLeptonDrSeparationFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genMatchedJetLeptonDrSeparation","genJetsPassingDrSeparationCut"),
		minNumber = cms.uint32(2)
		)

genMatchedLeadingEleDrSeparationFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genMatchedJetLeptonDrSeparation","matchedLeadingGenElePassingDrSeparationCut"),
		minNumber = cms.uint32(1)
		)

genMatchedSubleadingEleDrSeparationFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genMatchedJetLeptonDrSeparation","matchedSubleadingGenElePassingDrSeparationCut"),
		minNumber = cms.uint32(1)
		)

genMatchedJetLeptonDrSeparationSeq = cms.Sequence(
		genMatchedJetLeptonDrSeparation
		*genMatchedJetLeptonDrSeparationFilter
		*genMatchedLeadingEleDrSeparationFilter
		*genMatchedSubleadingEleDrSeparationFilter
		)

## end modules which apply dR(lepton, jet) cut


## modules which slim the leptons from dR(l,j) filter down to a single module name, instead of collections identified with two strings 
pickGenMatchedLeadingEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("genMatchedJetLeptonDrSeparation","matchedLeadingGenElePassingDrSeparationCut"),
		cut = cms.string("")
		)

pickGenMatchedSubleadingEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("genMatchedJetLeptonDrSeparation","matchedSubleadingGenElePassingDrSeparationCut"),
		cut = cms.string("")
		)

pickGenMatchedEleSeq = cms.Sequence(pickGenMatchedLeadingEle*pickGenMatchedSubleadingEle)

## end modules which slim leptons from dR(l,j) filter down to a signle module name


## modules which require one gen electron have pt>60
mergeGenMatchedEles = cms.EDProducer("CandViewMerger",
		src = cms.VInputTag("pickGenMatchedLeadingEle","pickGenMatchedSubleadingEle")
		)

requireGenMatchedHighPtEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("mergeGenMatchedEles"),
		cut = cms.string("pt>60")
		)

requireGenMatchedHighPtEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("requireGenMatchedHighPtEle"),
		minNumber = cms.uint32(1)
		)

requireGenMatchedHighPtEleSeq = cms.Sequence(
		mergeGenMatchedEles
		*requireGenMatchedHighPtEle
		*requireGenMatchedHighPtEleFilter
		)

## end modules which require one gen electron have pt>60

## modules which apply the dilepton mass cut to electrons which are being studied with gen jets
genMatchedDiLeptonCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("pickGenMatchedLeadingEle pickGenMatchedSubleadingEle"),
		role = cms.string("matchedLeadingEle matchedSubleadingEle"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 200")
		)

genMatchedDiLeptonCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genMatchedDiLeptonCandidate"),
		minNumber = cms.uint32(1)
		)

genMatchedDiLeptonCandidateSeq = cms.Sequence(genMatchedDiLeptonCandidate*genMatchedDiLeptonCandidateFilter)

## end modules which apply the dilepton mass cut


## modules which apply the four object mass cut using gen electrons and gen jets matched to WR decay products
genMatchedFourObjMass = cms.EDProducer("applyFourObjMassCutTwoInputLeptonColls",
		outputJetsCollectionName = cms.string("genJetsPassingFourObjMassCut"),
		outputLeptonsOneCollectionName = cms.string("genMatchedLeadingElectronPassingFourObjMassCut"),
		outputLeptonsTwoCollectionName = cms.string("genMatchedSubleadingElectronPassingFourObjMassCut"),
		minFourObjMassCut = cms.double(600.0),
		minDileptonMassCut = cms.double(200.0),
		minLeptonDrSeparation = cms.double(0.4),  #separation btwn selected leptons
		inputJetsCollTag = cms.InputTag("genMatchedJetLeptonDrSeparation","genJetsPassingDrSeparationCut"),
		inputLeptonsOneCollTag = cms.InputTag("pickGenMatchedLeadingEle"),
		inputLeptonsTwoCollTag = cms.InputTag("pickGenMatchedSubleadingEle"),
		)

genMatchedFourObjMassLeptonsOneFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genMatchedFourObjMass","genMatchedLeadingElectronPassingFourObjMassCut"),
		minNumber = cms.uint32(1)
		)

genMatchedFourObjMassLeptonsTwoFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genMatchedFourObjMass","genMatchedSubleadingElectronPassingFourObjMassCut"),
		minNumber = cms.uint32(1)
		)

genMatchedFourObjMassJetsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genMatchedFourObjMass","genJetsPassingFourObjMassCut"),
		minNumber = cms.uint32(2)
		)

genMatchedFourObjMassSeq = cms.Sequence(
		genMatchedFourObjMass
		*genMatchedFourObjMassLeptonsOneFilter
		*genMatchedFourObjMassLeptonsTwoFilter
		*genMatchedFourObjMassJetsFilter
		)

## end modules which apply the four object mass cut using gen electrons and gen jets matched to WR decay products


##filters on pt and eta of gen leptons and jets
ptEtaRestrictedGenJet = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("slimmedGenJets"),
		src = cms.InputTag("bareGenJet"),
		cut = cms.string("pt>40 && abs(eta) < 2.5")
		)

ptEtaRestrictedGenJetFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedGenJet"),
		minNumber = cms.uint32(2)
		)

ptEtaRestrictedSubleadingGenEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareGenEle"),
		cut = cms.string("pt>40 && abs(eta) < 2.5")
		)

ptEtaRestrictedSubleadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedSubleadingGenEle"),
		minNumber = cms.uint32(2)
		)

ptEtaRestrictedLeadingGenEle = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("ptEtaRestrictedSubleadingGenEle"),
		cut = cms.string("pt>60")
		)

ptEtaRestrictedLeadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedLeadingGenEle"),
		minNumber = cms.uint32(1)
		)

ptEtaConstrainedSeq = cms.Sequence(ptEtaRestrictedGenJet*ptEtaRestrictedGenJetFilter*ptEtaRestrictedSubleadingGenEle
		*ptEtaRestrictedSubleadingGenEleFilter*ptEtaRestrictedLeadingGenEle*ptEtaRestrictedLeadingGenEleFilter)


#these filters and sequence look for gen quarks and leptons which are matched (by pdgId) to particles in the WR
#decay chain

ptEtaRestrictedMatchedGenQuark = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("etaRestrictedMatchedGenQuark"),
		#src = cms.InputTag("bareMatchedGenQuark"),
		cut = cms.string("pt>40 && abs(eta) < 2.5")
		)

ptEtaRestrictedMatchedGenQuarkFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedMatchedGenQuark"),
		minNumber = cms.uint32(2)
		)

ptEtaRestrictedMatchedSubleadingGenEle = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("bareMatchedSubleadingGenEle"),
		src = cms.InputTag("etaRestrictedMatchedGenSubleadingLepton"),
		cut = cms.string("pt>40 && abs(eta) < 2.5")
		)

ptEtaRestrictedMatchedSubleadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedMatchedSubleadingGenEle"),
		minNumber = cms.uint32(1)
		)

ptEtaRestrictedMatchedLeadingGenEle = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("bareMatchedLeadingGenEle"),
		src = cms.InputTag("etaRestrictedMatchedGenLeadingLepton"),
		cut = cms.string("pt>60 && abs(eta) < 2.5")
		)

ptEtaRestrictedMatchedLeadingGenEleFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedMatchedLeadingGenEle"),
		minNumber = cms.uint32(1)
		)

matchedPtEtaConstrainedSeq = cms.Sequence(ptEtaRestrictedMatchedGenQuark*ptEtaRestrictedMatchedGenQuarkFilter
		*ptEtaRestrictedMatchedSubleadingGenEle*ptEtaRestrictedMatchedSubleadingGenEleFilter
		*ptEtaRestrictedMatchedLeadingGenEle*ptEtaRestrictedMatchedLeadingGenEleFilter)



##filters on dilepton mass, and the producer necessary to facilitate this filter
genDiElectronCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("ptEtaRestrictedLeadingGenEle ptEtaRestrictedSubleadingGenEle"),
		role = cms.string("leading subleading"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 200 && daughter(0).pt > daughter(1).pt")
		)

genDiElectronCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genDiElectronCandidate"),
		minNumber = cms.uint32(1)
		)

genDiElectronCandidateSeq = cms.Sequence(genDiElectronCandidate*genDiElectronCandidateFilter)

genMatchedDiElectronCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("ptEtaRestrictedMatchedLeadingGenEle ptEtaRestrictedMatchedSubleadingGenEle"),
		role = cms.string("matchedLeadingEle matchedSubleadingEle"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 200 && daughter(0).pt > daughter(1).pt")
		)

genMatchedDiElectronCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("genMatchedDiElectronCandidate"),
		minNumber = cms.uint32(1)
		)

genMatchedDiElectronCandidateSeq = cms.Sequence(genMatchedDiElectronCandidate*genMatchedDiElectronCandidateFilter)



