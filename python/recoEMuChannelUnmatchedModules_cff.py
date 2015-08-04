import FWCore.ParameterSet.Config as cms

#this is a copy of recoElectronChannelUnmatchedModules_cff

## transform the objects in slimmedJets, slimmedElectrons, and slimmedMuons into reco::Candidate objects
emuBareRecoJet = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedJets"),
		cut = cms.string("abs(eta) < 2.5 && (neutralHadronEnergyFraction<0.90 && neutralEmEnergyFraction<0.9 && (chargedMultiplicity+neutralMultiplicity)>1 && muonEnergyFraction<0.8) && ((abs(eta)<=2.4 && chargedHadronEnergyFraction>0 && chargedMultiplicity>0 && chargedEmEnergyFraction<0.90) || abs(eta)>2.4) ")
		)

emuBareRecoJetFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuBareRecoJet"),
		minNumber = cms.uint32(2)
		)

emuBareRecoLeptonOne = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedElectrons"),
		cut = cms.string("abs(eta) < 2.5 && pt>20")
		)

emuBareRecoLeptonOneFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuBareRecoLeptonOne"),
		minNumber = cms.uint32(1)
		)

emuBareRecoLeptonTwo = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedMuons"),
		cut = cms.string("abs(eta) < 2.4 && pt>20")
		)

emuBareRecoLeptonTwoFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuBareRecoLeptonTwo"),
		minNumber = cms.uint32(1)
		)

emuMergeLeptons = cms.EDProducer("CandViewMerger",
		src = cms.VInputTag("emuBareRecoLeptonOne","emuBareRecoLeptonTwo")
		)

emuRequireHighPtLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("emuMergeLeptons"),
		cut = cms.string("pt>40")
		)

emuRequireHighPtLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuRequireHighPtLepton"),
		minNumber = cms.uint32(1)
		)

emuBareRecoParticleSeq = cms.Sequence(
		emuBareRecoJet
		*emuBareRecoJetFilter
		*emuBareRecoLeptonOne
		*emuBareRecoLeptonOneFilter
		*emuBareRecoLeptonTwo
		*emuBareRecoLeptonTwoFilter
		*emuMergeLeptons
		*emuRequireHighPtLepton
		*emuRequireHighPtLeptonFilter
		)

##########################
##these modules apply the dR separation cut btwn leptons and jets after they
##are selected by the bareReco modules

emuBareRecoJetLeptonDrSeparation = cms.EDProducer("applyLeptonJetDrCutMixedLeptonFlavor",
		outputJetsCollectionName = cms.string("bareJetsPassingDrSeparationCut"),
		minDrSeparation = cms.double(0.4),
		inputJetsCollTag = cms.InputTag("emuBareRecoJet"),
		inputLeptonsOneCollTag = cms.InputTag("emuBareRecoLeptonOne"),
		inputLeptonsTwoCollTag = cms.InputTag("emuBareRecoLeptonTwo"),
		minDileptonMassCut = cms.double(-1),
		)

emuBareRecoJetLeptonDrSeparationFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuBareRecoJetLeptonDrSeparation","bareJetsPassingDrSeparationCut"),
		minNumber = cms.uint32(2)
		)

emuBareRecoDrSeparationSeq = cms.Sequence(
		emuBareRecoJetLeptonDrSeparation
		*emuBareRecoJetLeptonDrSeparationFilter
		)
##########################



