import FWCore.ParameterSet.Config as cms

## checking Zee events
## run this after triggerFilter and HEEPID selector
zeeCheckLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("HEEPIDSelector"),
		#src = cms.InputTag("slimmedElectrons"),
		cut = cms.string("")
		)

zeeCheckLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("zeeCheckLepton"),
		minNumber = cms.uint32(2)
		)

checkZeeSeq = cms.Sequence(
		zeeCheckLepton
		*zeeCheckLeptonFilter
		)

## transform the objects in slimmedJets and slimmedElectrons into reco::Candidate objects
bareRecoJet = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("slimmedJets"),
		cut = cms.string("pt>40 && abs(eta) < 2.5 && (neutralHadronEnergyFraction<0.90 && neutralEmEnergyFraction<0.9 && (chargedMultiplicity+neutralMultiplicity)>1 && muonEnergyFraction<0.8) && ((abs(eta)<=2.4 && chargedHadronEnergyFraction>0 && chargedMultiplicity>0 && chargedEmEnergyFraction<0.90) || abs(eta)>2.4) ")
		)

bareRecoJetFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareRecoJet"),
		minNumber = cms.uint32(2)
		)

bareRecoLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("HEEPIDSelector"),
		cut = cms.string("pt > 40")
		)

bareRecoLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareRecoLepton"),
		minNumber = cms.uint32(2)
		)

bareRecoParticleSeq = cms.Sequence(bareRecoJet*bareRecoJetFilter*bareRecoLepton*bareRecoLeptonFilter)

## for the low mass sideband, skip evts with two HEEP electrons with dilepton mass > 200 GeV

eeEarlyRecoLowWrMassFilter = cms.EDFilter("hasNoHighMassWrObjects",
		maxWrMass = cms.double(14000.0),
		maxDileptonMass = cms.double(200.0),
		inputLeadLeptonsCollTag = cms.InputTag("bareRecoLepton"),
		inputSubleadLeptonsCollTag = cms.InputTag("bareRecoLepton"),
		inputJetsCollTag = cms.InputTag("bareRecoJet")
		)

eeEarlyLowDileptonMassSeq = cms.Sequence(eeEarlyRecoLowWrMassFilter)

##########################
##these modules apply the dR separation cut btwn leptons and jets just after they
##are selected by the bareReco modules

bareRecoJetLeptonDrSeparation = cms.EDProducer("applyLeptonJetDrCut",
		outputLeptonsCollectionName = cms.string("bareLeptonsPassingDrSeparationCut"),
		minDrSeparation = cms.double(0.4),
		inputJetsCollTag = cms.InputTag("bareRecoJet"),
		inputLeptonsCollTag = cms.InputTag("bareRecoLepton"),
		minDileptonMassCut = cms.double(-1),
		)

bareRecoJetLeptonDrSeparationFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("bareRecoJetLeptonDrSeparation","bareLeptonsPassingDrSeparationCut"),
		minNumber = cms.uint32(2)
		)

bareRecoDrSeparationSeq = cms.Sequence(bareRecoJetLeptonDrSeparation*bareRecoJetLeptonDrSeparationFilter)
##########################


## apply pT and eta cuts to the bare reco::Candidate leptons and jets
ptEtaRestrictedRecoLeptons = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("bareRecoLepton"),
		cut = cms.string("abs(eta) < 2.5 && pt>40")
		)
ptEtaRestrictedRecoLeptonsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		minNumber = cms.uint32(2)
		)

ptEtaRestrictedLeadRecoLepton = cms.EDFilter("CandViewSelector",
		src = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		cut = cms.string("pt>60")
		)
ptEtaRestrictedLeadRecoLeptonFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedLeadRecoLepton"),
		minNumber = cms.uint32(1)
		)

ptEtaRestrictedRecoJets = cms.EDFilter("CandViewSelector",
		#src = cms.InputTag("bareRecoJetLeptonDrSeparation","bareJetsPassingDrSeparationCut"),
		src = cms.InputTag("bareRecoJet"),
		cut = cms.string("abs(eta) < 2.5 && pt>40")
		)
ptEtaRestrictedRecoJetsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("ptEtaRestrictedRecoJets"),
		minNumber = cms.uint32(2)
		)

ptEtaRestrictedSeq = cms.Sequence(
		ptEtaRestrictedRecoLeptons
		*ptEtaRestrictedRecoLeptonsFilter
		*ptEtaRestrictedLeadRecoLepton
		*ptEtaRestrictedLeadRecoLeptonFilter
		*ptEtaRestrictedRecoJets
		*ptEtaRestrictedRecoJetsFilter
		)
## end modules which apply pT and eta cuts to lepton and jets

## require that no LLJJ object appear with mass > 600 GeV in the evt
recoDiLeptonCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("ptEtaRestrictedLeadRecoLepton ptEtaRestrictedRecoLeptons"),
		role = cms.string("leadingLepton subleadingLepton"),
		checkCharge = cms.bool(False),
		cut = cms.string("0 < mass < 200 && daughter(0).pt > daughter(1).pt")
		)

recoDiLeptonCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("recoDiLeptonCandidate"),
		minNumber = cms.uint32(1)
		)

recoLowWrMassFilter = cms.EDFilter("hasNoHighMassWrObjects",
		maxWrMass = cms.double(600.0),
		maxDileptonMass = cms.double(200.0),
		inputLeadLeptonsCollTag = cms.InputTag("ptEtaRestrictedLeadRecoLepton"),
		inputSubleadLeptonsCollTag = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		inputJetsCollTag = cms.InputTag("ptEtaRestrictedRecoJets")
		)

lowMassLLJJObjectSeq = cms.Sequence(
		recoDiLeptonCandidate
		*recoDiLeptonCandidateFilter
		*recoLowWrMassFilter
		)
## end low mass LLJJ filter sequence


## these modules apply the dilepton mass cut to matched reco leptons
recoDileptonCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("ptEtaRestrictedLeadRecoLepton ptEtaRestrictedRecoLeptons"),
		role = cms.string("leadingLepton subleadingLepton"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 200 && daughter(0).pt > daughter(1).pt")
		)
recoDileptonCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("recoDileptonCandidate"),
		minNumber = cms.uint32(1)
		)

recoDileptonCandidateSeq = cms.Sequence(recoDileptonCandidate*recoDileptonCandidateFilter)

## end modules which apply dilepton mass cut to reco leptons


## these modules apply the dR(lepton,jet) separation cuts to jets which have
## passed the pT and eta filters
recoDrSeparation = cms.EDProducer("applyLeptonJetDrCut",
		outputJetsCollectionName = cms.string("jetsPassingDrSeparationCut"),
		minDrSeparation = cms.double(0.4),
		inputJetsCollTag = cms.InputTag("ptEtaRestrictedRecoJets"),
		inputLeptonsCollTag = cms.InputTag("ptEtaRestrictedRecoLeptons"),
		minDileptonMassCut = cms.double(200)
		)

recoDrSeparationFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("recoDrSeparation","jetsPassingDrSeparationCut"),
		minNumber = cms.uint32(2)
		)

recoDrSeparationSeq = cms.Sequence(recoDrSeparation*recoDrSeparationFilter)
## end modules which apply dR(lepton,jet) separation cuts


## these modules apply the four object inv mass cut
## if the two leading leptons and jets satisfy the four obj inv mass cut, then the
## producer adds two new collections - one for the two leptons, the other for the
## two jets - to the event
## if these two collections exist, and two objects exist in each collection, then
## the four obj mass cut is satisfied
recoFourObjMass = cms.EDProducer("applyFourObjMassCut",
		outputJetsCollectionName = cms.string("leadingJetsPassingFourObjMassCut"),
		outputLeptonsCollectionName = cms.string("leadingLeptonsPassingFourObjMassCut"),
		minFourObjMassCut = cms.double(600.0),
		minDileptonMassCut = cms.double(200.0),
		inputJetsCollTag = cms.InputTag("recoDrSeparation","jetsPassingDrSeparationCut"),
		inputLeptonsCollTag = cms.InputTag("ptEtaRestrictedRecoLeptons")
		)

recoFourObjMassLeptonsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("recoFourObjMass","leadingLeptonsPassingFourObjMassCut"),
		minNumber = cms.uint32(2)
		)

recoFourObjMassJetsFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("recoFourObjMass","leadingJetsPassingFourObjMassCut"),
		minNumber = cms.uint32(2)
		)

recoFourObjMassSeq = cms.Sequence(
		recoFourObjMass
		*recoFourObjMassLeptonsFilter
		*recoFourObjMassJetsFilter
		)

## end modules which apply the four object invariant mass cut


