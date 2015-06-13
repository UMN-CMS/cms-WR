import FWCore.ParameterSet.Config as cms

# Selection of electrons from the slimmedElectrons collection (PAT)
""" \addtogroup electronSkim_Group electronSkim sequences
@{
"""

### make sure the evt has at least two jets, and one has a nontrivial pT
wRhardJet = cms.EDFilter("PATJetRefSelector",
		src = cms.InputTag("slimmedJets"),
		cut = cms.string("pt>20")
		)

wRsoftJet = cms.EDFilter("PATJetRefSelector",
		src = cms.InputTag("slimmedJets"),
		cut = cms.string("pt>8")
		)

wRdiJetCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("wRhardJet wRsoftJet"),
		role = cms.string("leadingJet subleadingJet"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 0")
		)

wRdiJetCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("wRdiJetCandidate"),
		minNumber = cms.uint32(1)
		)

### select leading electron \ingroup electronSkim_Group
wRleadingElectron = cms.EDFilter("PATElectronRefSelector",
                                 src = cms.InputTag("slimmedElectrons"),
                                 cut = cms.string("pt>40"),
                                 )

### select subleading electron
wRsubleadingElectron = cms.EDFilter("PATElectronRefSelector",
                                 src = cms.InputTag("slimmedElectrons"),
                                 cut = cms.string("pt>20"),
                                 )
#wRpreSelectedElectrons = cms.EDProducer("CandViewMerger",

### create di-electron pair in signal region
wRdiElectronCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                       decay = cms.string("wRleadingElectron wRsubleadingElectron"),
                                       role = cms.string("leading subleading"),
                                       checkCharge = cms.bool(False),
                                       # the cut on the pt of the daughter is to respect the order of leading and subleading:
                                           # if both electrons have pt>60 GeV there will be two di-electron candidates with inversed order
                                       #cut = cms.string("mass > 200 && daughter(0).pt>daughter(1).pt"),
                                       cut = cms.string("mass > 0"),
									   )

### filter: at least one di-electron candidate in signal region
wRdiElectronCandidateFilter = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("wRdiElectronCandidate"),
                                           minNumber = cms.uint32(1)
                                           )

### create an object from two jets and two electrons in the evt, and cut on its mass
wRdiLeptonDijetCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
		decay = cms.string("wRdiJetCandidate wRdiElectronCandidate"),
		role = cms.string("dijet dilepton"),
		checkCharge = cms.bool(False),
		cut = cms.string("mass > 500")
		)

### filter: require at least one LLJJ object in the evt
wRdiLeptonDijetCandidateFilter = cms.EDFilter("CandViewCountFilter",
		src = cms.InputTag("wRdiLeptonDijetCandidate"),
		minNumber = cms.uint32(1)
		)


### create di-electron pair in sideband region
wRdiElectronSidebandCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                       decay = cms.string("wRleadingElectron wRsubleadingElectron"),
                                       checkCharge = cms.bool(False),
                                       cut = cms.string("0< mass < 200 && daughter(0).pt>daughter(1).pt")
                                       )

### filter: at least one di-electron candidate in sideband region
wRdiElectronSidebandCandidateFilter = cms.EDFilter("CandViewCountFilter",
                                                   src = cms.InputTag("wRdiElectronSidebandCandidate"),
                                                   minNumber = cms.uint32(1)
                                                   )

### di-jet selection sequence
wRjetSelectionSeq = cms.Sequence(wRhardJet + wRsoftJet)
### di-electron selection sequence
wRelectronSelectionSeq = cms.Sequence(wRleadingElectron + wRsubleadingElectron)
### di-electron selection in signal region sequence
wRdiElectronSignalSeq = cms.Sequence(wRelectronSelectionSeq * wRdiElectronCandidate * wRdiElectronCandidateFilter)
### di-electron and four object selection in signal region sequence
wRdiElectronAndFourObjSignalSeq = cms.Sequence(
		wRjetSelectionSeq
		*wRdiElectronSignalSeq
		*wRdiJetCandidate
		*wRdiJetCandidateFilter
		*wRdiLeptonDijetCandidate
		*wRdiLeptonDijetCandidateFilter
		)
### di-electron selection in sideband region sequence
wRdiElectronSidebandSeq = cms.Sequence(wRelectronSelectionSeq * wRdiElectronSidebandCandidate * wRdiElectronSidebandCandidateFilter)


### @}
