import FWCore.ParameterSet.Config as cms

# Selection of muons from the slimmedMuons collection (PAT)
""" \addtogroup muonSkim_Group muonSkim sequences
@{
"""

### make sure the evt has at least two jets, and one has a nontrivial pT
muwRhardJet = cms.EDFilter("PATJetRefSelector",
		src = cms.InputTag("slimmedJets"),
		cut = cms.string("pt>20")
		)

muwRsoftJet = cms.EDFilter("PATJetRefSelector",
		src = cms.InputTag("slimmedJets"),
		cut = cms.string("pt>8")
		)

muwRdiJetCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("muwRhardJet muwRsoftJet"),
		role = cms.string("leadingJet subleadingJet"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 0")
		)

muwRdiJetCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("muwRdiJetCandidate"),
		minNumber = cms.uint32(1)
		)



### select leading muon \ingroup muonSkim_Group
wRleadingMuon = cms.EDFilter("PATMuonRefSelector",
                             src = cms.InputTag("slimmedMuons"),
                             cut = cms.string("pt>45"),
                         )

### select subleading muon
wRsubleadingMuon = cms.EDFilter("PATMuonRefSelector",
                                src = cms.InputTag("slimmedMuons"),
                                cut = cms.string("pt>30"),
                            )
### select loose-ID jets
#wRlooseJet = cms.EDFilter("PATJetSelector",
#                            src = cms.InputTag("slimmedJets"),
#                            cut = cms.string("(neutralHadronEnergyFraction<0.90 && neutralEmEnergyFraction<0.9 && (chargedMultiplicity+neutralMultiplicity)>1 && muonEnergyFraction<0.8) && ((abs(eta)<=2.4 && chargedHadronEnergyFraction>0 && chargedMultiplicity>0 && chargedEmEnergyFraction<0.90) || abs(eta)>2.4)"),
#                            #cut = cms.string("pt>60"),
#                            )

### create di-muon pair in signal region
wRdiMuonCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string("wRleadingMuon wRsubleadingMuon"),
                                   role = cms.string("leading subleading"),
                                   checkCharge = cms.bool(False),
                                   # the cut on the pt of the daughter is to respect the order of leading and subleading:
                                   # if both muons have pt>60 GeV there will be two di-muon candidates with inversed order
                                   #cut = cms.string("mass > 200 && daughter(0).pt>daughter(1).pt"),
                                   cut = cms.string("mass > 0  && daughter(0).pt>daughter(1).pt"),
                                   
									   )

### filter: at least one di-muon candidate in signal region
wRdiMuonCandidateFilter = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("wRdiMuonCandidate"),
                                           minNumber = cms.uint32(1)
                                           )


### create an object from two jets and two muons in the evt, and cut on its mass
wRdiMuonDijetCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("muwRdiJetCandidate wRdiMuonCandidate"),
		role = cms.string("dijet dilepton"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 500")
		)

### filter: require at least one LLJJ object in the evt
wRdiMuonDijetCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("wRdiMuonDijetCandidate"),
		minNumber = cms.uint32(1)
		)



### create di-muon pair in sideband region
wRdiMuonSidebandCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                       decay = cms.string("wRleadingMuon wRsubleadingMuon"),
                                       checkCharge = cms.bool(False),
                                       cut = cms.string("0< mass < 200 && daughter(0).pt>daughter(1).pt")
                                       )

### filter: at least one di-muon candidate in sideband region
wRdiMuonSidebandCandidateFilter = cms.EDFilter("CandViewCountFilter",
                                                   src = cms.InputTag("wRdiMuonSidebandCandidate"),
                                                   minNumber = cms.uint32(1)
                                                   )

### di-jet selection sequence
muwRjetSelectionSeq = cms.Sequence(muwRhardJet + muwRsoftJet)
### di-muon selection sequence
wRmuonSelectionSeq = cms.Sequence(wRleadingMuon + wRsubleadingMuon)
### di-muon selection in signal region sequence
wRdiMuonSignalSeq = cms.Sequence(wRmuonSelectionSeq * wRdiMuonCandidate * wRdiMuonCandidateFilter)
### di-muon and four object selection in signal region sequence
wRdiMuonAndFourObjSignalSeq = cms.Sequence(
		muwRjetSelectionSeq
		*wRdiMuonSignalSeq
		*muwRdiJetCandidate
		*muwRdiJetCandidateFilter
		*wRdiMuonDijetCandidate
		*wRdiMuonDijetCandidateFilter
		)
### di-muon selection in sideband region sequence
wRdiMuonSidebandSeq = cms.Sequence(wRmuonSelectionSeq * wRdiMuonSidebandCandidate * wRdiMuonSidebandCandidateFilter)

### di-muon and low four object mass selection sequence
wRdiMuonAndLowMassSignalSeq = cms.Sequence(
		muwRjetSelectionSeq
		*wRdiMuonSignalSeq
		*muwRdiJetCandidate
		*muwRdiJetCandidateFilter
		*wRdiMuonDijetCandidate
		*~wRdiMuonDijetCandidateFilter
		)

### @}
