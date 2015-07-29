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

emuBareRecoLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedElectrons"),
		cut = cms.string("")
		)

emuBareRecoLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuBareRecoLepton"),
		minNumber = cms.uint32(1)
		)

emuBareRecoLeptonTwo = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedMuons"),
		cut = cms.string("")
		)

emuBareRecoLeptonTwoFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuBareRecoLeptonTwo"),
		minNumber = cms.uint32(1)
		)

emuCombineBareRecoLeptonsIntoOneCollection = cms.EDProducer('combineLeptonsProducer',
		inputLeptonOneTag = cms.InputTag("emuBareRecoLepton"),
		inputLeptonTwoTag = cms.InputTag("emuBareRecoLeptonTwo"),
		outputLeptonCollName = cms.string("bareRecoElectronsAndMuons")
		)

emuBareRecoParticleSeq = cms.Sequence(
		emuBareRecoJet
		*emuBareRecoJetFilter
		*emuBareRecoLepton
		*emuBareRecoLeptonFilter
		*emuBareRecoLeptonTwo
		*emuBareRecoLeptonTwoFilter
		#*emuCombineBareRecoLeptonsIntoOneCollection
		)

