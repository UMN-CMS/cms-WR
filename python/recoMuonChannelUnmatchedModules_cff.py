import FWCore.ParameterSet.Config as cms

#this is a copy of recoElectronChannelUnmatchedModules_cff

## transform the objects in slimmedJets and slimmedMuons into reco::Candidate objects
mumuBareRecoJet = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedJets"),
		cut = cms.string("abs(eta) < 2.5 && (neutralHadronEnergyFraction<0.90 && neutralEmEnergyFraction<0.9 && (chargedMultiplicity+neutralMultiplicity)>1 && muonEnergyFraction<0.8) && ((abs(eta)<=2.4 && chargedHadronEnergyFraction>0 && chargedMultiplicity>0 && chargedEmEnergyFraction<0.90) || abs(eta)>2.4) ")
		)

mumuBareRecoJetFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("mumuBareRecoJet"),
		minNumber = cms.uint32(2)
		)

mumuBareRecoLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedMuons"),
		cut = cms.string("abs(eta) < 2.4 && pt>20")
		)

mumuBareRecoLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("mumuBareRecoLepton"),
		minNumber = cms.uint32(2)
		)

mumuBareRecoParticleSeq = cms.Sequence(mumuBareRecoJet*mumuBareRecoJetFilter*mumuBareRecoLepton*mumuBareRecoLeptonFilter)

##########################
##these modules apply the dR separation cut btwn leptons and jets after they
##are selected by the bareReco modules

mumuBareRecoJetLeptonDrSeparation = cms.EDProducer("applyLeptonJetDrCut",
		outputJetsCollectionName = cms.string("bareJetsPassingDrSeparationCut"),
		minDrSeparation = cms.double(0.4),
		inputJetsCollTag = cms.InputTag("mumuBareRecoJet"),
		inputLeptonsCollTag = cms.InputTag("mumuBareRecoLepton"),
		minDileptonMassCut = cms.double(-1),
		)

mumuBareRecoJetLeptonDrSeparationFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("mumuBareRecoJetLeptonDrSeparation","bareJetsPassingDrSeparationCut"),
		minNumber = cms.uint32(2)
		)

mumuBareRecoDrSeparationSeq = cms.Sequence(mumuBareRecoJetLeptonDrSeparation*mumuBareRecoJetLeptonDrSeparationFilter)
##########################


