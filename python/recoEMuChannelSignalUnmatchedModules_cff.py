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
		minPt = cms.double(50)
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
		src = cms.InputTag("HEEPIDSelector"),
		cut = cms.string("pt>40 && abs(eta) < 2.5")
		)

emuBareRecoLeptonOneFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuBareRecoLeptonOne"),
		minNumber = cms.uint32(1)
		)

emuBareRecoLeptonTwo = cms.EDFilter("CandViewSelector",
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

## for the low dilepton mass sideband, skip evts which have leptons passing ID requirements
## which also have dilepton mass > 200 GeV
#
#emuEarlyRecoLowMassFilter = cms.EDFilter("hasNoHighMassWrObjects",
#		maxWrMass = cms.double(14000.0),
#		maxDileptonMass = cms.double(200.0),
#		inputLeadLeptonsCollTag = cms.InputTag("emuBareRecoLeptonOne"),
#		inputSubleadLeptonsCollTag = cms.InputTag("emuBareRecoLeptonTwo"),
#		inputJetsCollTag = cms.InputTag("emuBareRecoJet")
#		)
#
#emuEarlyLowDileptonMassSeq = cms.Sequence(emuEarlyRecoLowMassFilter)

## end modules to skip high dilepton mass evts

##########################
##these modules apply the dR separation cut btwn leptons and jets after they
##are selected by the bareReco modules

emuBareRecoJetLeptonDrSeparation = cms.EDProducer("applyLeptonJetDrCutMixedLeptonFlavor",
		outputJetsCollectionName = cms.string("jetsPassingDrSeparationCut"),
		outputLeptonsOneCollectionName = cms.string("leptonsOnePassingDrSeparationCut"),
		outputLeptonsTwoCollectionName = cms.string("leptonsTwoPassingDrSeparationCut"),
		minDrSeparation = cms.double(0.4),
		inputJetsCollTag = cms.InputTag("emuBareRecoJet"),
		inputLeptonsOneCollTag = cms.InputTag("emuBareRecoLeptonOne"),
		inputLeptonsTwoCollTag = cms.InputTag("emuBareRecoLeptonTwo"),
		minDileptonMassCut = cms.double(-1),
		minLeptonDrSeparation = cms.double(0.4)
		)

emuBareRecoJetLeptonDrSeparationJetsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuBareRecoJetLeptonDrSeparation","jetsPassingDrSeparationCut"),
		minNumber = cms.uint32(2)
		)

emuBareRecoJetLeptonDrSeparationLeptonsOneFilter = emuBareRecoJetLeptonDrSeparationJetsFilter.clone(
		src = cms.InputTag("emuBareRecoJetLeptonDrSeparation","leptonsOnePassingDrSeparationCut"),
		minNumber = cms.uint32(1)
		)

emuBareRecoJetLeptonDrSeparationLeptonsTwoFilter = emuBareRecoJetLeptonDrSeparationLeptonsOneFilter.clone(
		src = cms.InputTag("emuBareRecoJetLeptonDrSeparation","leptonsTwoPassingDrSeparationCut")
		)

emuBareRecoDrSeparationSeq = cms.Sequence(
		emuBareRecoJetLeptonDrSeparation
		*emuBareRecoJetLeptonDrSeparationJetsFilter
		*emuBareRecoJetLeptonDrSeparationLeptonsOneFilter
		*emuBareRecoJetLeptonDrSeparationLeptonsTwoFilter
		)
#end modules which apply dR(L,J) cuts
##########################

### begin modules which collimate lepton collection names in preparation for high pT lepton cut
emuCollimateLeptonsOne = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("emuBareRecoJetLeptonDrSeparation","leptonsOnePassingDrSeparationCut"),
		cut = cms.string("")
		)

emuCollimateLeptonsTwo = emuCollimateLeptonsOne.clone(
		src = cms.InputTag("emuBareRecoJetLeptonDrSeparation","leptonsTwoPassingDrSeparationCut")
		)

emuCollimateLeptonsSeq = cms.Sequence(emuCollimateLeptonsOne*emuCollimateLeptonsTwo)

### end modules which collimate lepton collection names in preparation for high pT lepton cut


## the jets have already been required to pass the loose jet ID and have pt>40
## the leptons have already been required to have pt>40
## now require one lepton has pt>60
## run the low mass filter plugin in a module, but set the max lljj mass to 14000, and the max dilepton mass to 200

emuMergeLeptons = cms.EDProducer("CandViewMerger",
		src = cms.VInputTag("emuCollimateLeptonsOne","emuCollimateLeptonsTwo")
		)

emuRequireHighPtLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("emuMergeLeptons"),
		cut = cms.string("pt>60")
		)

emuRequireHighPtLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuRequireHighPtLepton"),
		minNumber = cms.uint32(1)
		)

#emuRecoLowMassFilter = cms.EDFilter("hasNoHighMassWrObjects",
#		maxWrMass = cms.double(14000.0),
#		maxDileptonMass = cms.double(200.0),
#		inputLeadLeptonsCollTag = cms.InputTag("emuBareRecoLeptonOne"),
#		inputSubleadLeptonsCollTag = cms.InputTag("emuBareRecoLeptonTwo"),
#		inputJetsCollTag = cms.InputTag("emuBareRecoJetLeptonDrSeparation","jetsPassingDrSeparationCut")
#		)
#
#emuLeadLeptonAndLowMassSeq = cms.Sequence(
#		emuMergeLeptons
#		*emuRequireHighPtLepton
#		*emuRequireHighPtLeptonFilter
#		*emuRecoLowMassFilter
#		)

## end lepton pt>60 and dilepton mass < 200 

## require one lepton has pT > 60
emuLeadLeptonSeq = cms.Sequence(
		emuMergeLeptons
		*emuRequireHighPtLepton
		*emuRequireHighPtLeptonFilter
		)
## end high pt lepton seq


## require the two leptons have dilepton mass > 200
emuRecoDiLeptonCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("emuCollimateLeptonsOne emuCollimateLeptonsTwo"),
		role = cms.string(""),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 200")
		)

emuRecoDiLeptonCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuRecoDiLeptonCandidate"),
		minNumber = cms.uint32(1)
		)

emuDileptonSignalSeq = cms.Sequence(
		emuRecoDiLeptonCandidate
		*emuRecoDiLeptonCandidateFilter
		)

##end dilepton mass requirement

##filter the jets and only retain those which are separated from the leptons
##by dR > 0.4
#emuFilteredRecoJetLeptonDrSeparation = cms.EDProducer("applyLeptonJetDrCutMixedLeptonFlavor",
#		outputJetsCollectionName = cms.string("filteredJetsPassingDrSeparationCut"),
#		minDrSeparation = cms.double(0.4),
#		inputJetsCollTag = cms.InputTag("emuBareRecoJet"),
#		inputLeptonsOneCollTag = cms.InputTag("emuBareRecoLeptonOne"),
#		inputLeptonsTwoCollTag = cms.InputTag("emuBareRecoLeptonTwo"),
#		minDileptonMassCut = cms.double(200),
#		minLeptonDrSeparation = cms.double(0.4)
#		)
#
#emuFilteredRecoJetLeptonDrSeparationFilter = cms.EDFilter("CandViewCountFilter",
#		src = cms.InputTag("emuFilteredRecoJetLeptonDrSeparation","filteredJetsPassingDrSeparationCut"),
#		minNumber = cms.uint32(2)
#		)
#
#emuFilteredRecoDrSeparationSeq = cms.Sequence(
#		emuFilteredRecoJetLeptonDrSeparation
#		*emuFilteredRecoJetLeptonDrSeparationFilter
#		)
#
##end dR(lepton,jet) filter modules


## apply four object mass cut
emuRecoFourObjMass = cms.EDProducer("applyFourObjMassCutTwoInputLeptonColls",
		outputJetsCollectionName = cms.string("jetsPassingFourObjMassCut"),
		outputLeptonsOneCollectionName = cms.string("leptonsOnePassingFourObjMassCut"),
		outputLeptonsTwoCollectionName = cms.string("leptonsTwoPassingFourObjMassCut"),
		minFourObjMassCut = cms.double(600.0),
		minDileptonMassCut = cms.double(200.0),
		minLeptonDrSeparation = cms.double(0.4),
		inputJetsCollTag = cms.InputTag("emuBareRecoJetLeptonDrSeparation","jetsPassingDrSeparationCut"),
		inputLeptonsOneCollTag = cms.InputTag("emuCollimateLeptonsOne"),
		inputLeptonsTwoCollTag = cms.InputTag("emuCollimateLeptonsTwo"),
		)

emuRecoFourObjMassLeptonsOneFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuRecoFourObjMass","leptonsOnePassingFourObjMassCut"),
		minNumber = cms.uint32(1)
		)

emuRecoFourObjMassLeptonsTwoFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuRecoFourObjMass","leptonsTwoPassingFourObjMassCut"),
		minNumber = cms.uint32(1)
		)

emuRecoFourObjMassJetsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("emuRecoFourObjMass","jetsPassingFourObjMassCut"),
		minNumber = cms.uint32(2)
		)

emuRecoFourObjMassSeq = cms.Sequence(
		emuRecoFourObjMass
		*emuRecoFourObjMassLeptonsOneFilter
		*emuRecoFourObjMassLeptonsTwoFilter
		*emuRecoFourObjMassJetsFilter
		)

## end modules which apply four object mass cut

