import FWCore.ParameterSet.Config as cms

# make a collection of TuneP muons which pass isHighPt ID
wrTunePMuProd = cms.EDProducer("TunePMuonProducer",
		src = cms.InputTag("slimmedMuons")
		)

wrTunePMuProdFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("wrTunePMuProd"),
		minNumber = cms.uint32(1)
		)

wrTunePMuProdSeq = cms.Sequence(
		wrTunePMuProd
		*wrTunePMuProdFilter
		)

# make a collection of TuneP muons which pass isHighPt ID
isHighPtMuProd = cms.EDProducer("produceIsHighPtMuons",
		src = cms.InputTag("wrTunePMuProd"),
		outputCollectionName = cms.string("TunePMuonsPassingIsHighPtID"),
		minPt = cms.double(45)
		)

isHighPtMuProdFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("isHighPtMuProd","TunePMuonsPassingIsHighPtID"),
		minNumber = cms.uint32(1)
		)

isHighPtMuSeq = cms.Sequence(isHighPtMuProd*isHighPtMuProdFilter)


#turn the TuneP muon and HEEP electron objects into reco::Candidate objects
#use this sequence to check the emu mass distribution, and ele and mu kinematics
electronCheck = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("HEEPIDSelector"),
		cut = cms.string("")
		)
electronCheckFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("electronCheck"),
		minNumber = cms.uint32(1)
		)

muonCheck = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("isHighPtMuProd","TunePMuonsPassingIsHighPtID"),
		cut = cms.string("abs(eta) < 2.4")
		)
muonCheckFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("muonCheck"),
		minNumber = cms.uint32(1)
		)

checkEMuSeq = cms.Sequence(
		electronCheck
		*electronCheckFilter
		*muonCheck
		*muonCheckFilter
		)


## transform the objects in slimmedJets, slimmedElectrons, and wrTunePMuProd into reco::Candidate objects
emuBareRecoJet = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedJets"),
		cut = cms.string("pt>40 && abs(eta) < 2.5 && (neutralHadronEnergyFraction<0.90 && neutralEmEnergyFraction<0.9 && (chargedMultiplicity+neutralMultiplicity)>1 && muonEnergyFraction<0.8) && ((abs(eta)<=2.4 && chargedHadronEnergyFraction>0 && chargedMultiplicity>0 && chargedEmEnergyFraction<0.90) || abs(eta)>2.4) ")
		)

emuBareRecoJetFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuBareRecoJet"),
		minNumber = cms.uint32(2)
		)

emuBareRecoLeptonOne = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("slimmedElectrons"),
		src = cms.InputTag("HEEPIDSelector"),
		#cut = cms.string("abs(eta) < 2.5 && pt>40")
		cut = cms.string("pt>40")
		)

emuBareRecoLeptonOneFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuBareRecoLeptonOne"),
		minNumber = cms.uint32(1)
		)

emuBareRecoLeptonTwo = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("wrTunePMuProd"),
		src = cms.InputTag("isHighPtMuProd","TunePMuonsPassingIsHighPtID"),
		cut = cms.string("abs(eta) < 2.4 && pt>50")  #pt>50 to get above isHighPt ID turnon curve
		)

emuBareRecoLeptonTwoFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuBareRecoLeptonTwo"),
		minNumber = cms.uint32(1)
		)

emuBareRecoParticleSeq = cms.Sequence(
		emuBareRecoJet
		*emuBareRecoJetFilter
		*emuBareRecoLeptonOne
		*emuBareRecoLeptonOneFilter
		*emuBareRecoLeptonTwo
		*emuBareRecoLeptonTwoFilter
		)

##########################
##these modules apply the dR separation cut btwn leptons and jets after they
##are selected by the bareReco modules

emuBareRecoJetLeptonDrSeparation = cms.EDProducer("applyLeptonJetDrCutMixedLeptonFlavor",
		outputJetsCollectionName = cms.string("jetsPassingDrSeparationCut"),
		minDrSeparation = cms.double(0.4),
		inputJetsCollTag = cms.InputTag("emuBareRecoJet"),
		inputLeptonsOneCollTag = cms.InputTag("emuBareRecoLeptonOne"),
		inputLeptonsTwoCollTag = cms.InputTag("emuBareRecoLeptonTwo"),
		minDileptonMassCut = cms.double(-1),
		)

emuBareRecoJetLeptonDrSeparationFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuBareRecoJetLeptonDrSeparation","jetsPassingDrSeparationCut"),
		minNumber = cms.uint32(2)
		)

emuBareRecoDrSeparationSeq = cms.Sequence(
		emuBareRecoJetLeptonDrSeparation
		*emuBareRecoJetLeptonDrSeparationFilter
		)
##########################

## the jets have already been required to pass the loose jet ID and have pt>40
## the leptons have already been required to have pt>40
## now require one lepton has pt>60
## run the low mass filter plugin in a module, but set the max lljj mass to 14000, and the max dilepton mass to 200

emuMergeLeptons = cms.EDProducer("CandViewMerger",
		src = cms.VInputTag("emuBareRecoLeptonOne","emuBareRecoLeptonTwo")
		)

emuRequireHighPtLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("emuMergeLeptons"),
		cut = cms.string("pt>60")
		)

emuRequireHighPtLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuRequireHighPtLepton"),
		minNumber = cms.uint32(1)
		)

emuRecoLowMassFilter = cms.EDFilter("hasNoHighMassWrObjects",
		maxWrMass = cms.double(14000.0),
		maxDileptonMass = cms.double(200.0),
		inputLeadLeptonsCollTag = cms.InputTag("emuBareRecoLeptonOne"),
		inputSubleadLeptonsCollTag = cms.InputTag("emuBareRecoLeptonTwo"),
		inputJetsCollTag = cms.InputTag("emuBareRecoJetLeptonDrSeparation","jetsPassingDrSeparationCut")
		)

emuLeadLeptonAndLowMassSeq = cms.Sequence(
		emuMergeLeptons
		*emuRequireHighPtLepton
		*emuRequireHighPtLeptonFilter
		*emuRecoLowMassFilter
		)

## end lepton pt>60 and dilepton mass < 200 


