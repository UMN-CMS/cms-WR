import FWCore.ParameterSet.Config as cms

# Selection of electrons from the slimmedElectrons collection (PAT)
""" \addtogroup electronSkim_Group electronSkim sequences
@{
"""


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
                                       cut = cms.string("mass > 200 && daughter(0).pt>daughter(1).pt"),
                                       )

### filter: at least one di-electron candidate in signal region
wRdiElectronCandidateFilter = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("wRdiElectronCandidate"),
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


### di-electron selection sequence
wRelectronSelectionSeq = cms.Sequence(wRleadingElectron + wRsubleadingElectron)
### di-electron selection in signal region sequence
wRdiElectronSignalSeq = cms.Sequence(wRelectronSelectionSeq * wRdiElectronCandidate * wRdiElectronCandidateFilter)
### di-electron selection in sideband region sequence
wRdiElectronSidebandSeq = cms.Sequence(wRelectronSelectionSeq * wRdiElectronSidebandCandidate * wRdiElectronSidebandCandidateFilter)


### @}
