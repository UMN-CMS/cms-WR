import FWCore.ParameterSet.Config as cms


############################################################
###### di electron selection
wRleadingElectron = cms.EDFilter("CandViewSelector",
                                 src = cms.InputTag("slimmedElectrons"),
                                 cut = cms.string("pt>60"),
                                 )

wRsubleadingElectron = cms.EDFilter("CandViewSelector",
                                 src = cms.InputTag("slimmedElectrons"),
                                 cut = cms.string("pt>40"),
                                 )
####### di electron signal
wRdiElectronCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                       decay = cms.string("wRleadingElectron wRsubleadingElectron"),
                                       role = cms.string("leading subleading"),
                                       checkCharge = cms.bool(False),
                                       # the cut on the pt of the daughter is to respect the order of leading and subleading:
                                           # if both electrons have pt>60 GeV there will be two di-electron candidates with inversed order
                                       cut = cms.string("mass > 200 && daughter(0).pt>daughter(1).pt"),
                                       )


wRdiElectronCandidateFilter = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("wRdiElectronCandidate"),
                                           minNumber = cms.uint32(1)
                                           )

####### di electron sideband
wRdiElectronSidebandCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                       decay = cms.string("wRleadingElectron wRsubleadingElectron"),
                                       checkCharge = cms.bool(False),
                                       cut = cms.string("0< mass < 200 && daughter(0).pt>daughter(1).pt")
                                       )

wRdiElectronSidebandCandidateFilter = cms.EDFilter("CandViewCountFilter",
                                                   src = cms.InputTag("wRdiElectronSidebandCandidate"),
                                                   minNumber = cms.uint32(1)
                                                   )


######## Sequences
wRelectronSelectionSeq = cms.Sequence(wRleadingElectron + wRsubleadingElectron)
wRdiElectronSignalSeq = cms.Sequence(wRelectronSelectionSeq * wRdiElectronCandidate * wRdiElectronCandidateFilter)
wRdiElectronSidebandSeq = cms.Sequence(wRelectronSelectionSeq * wRdiElectronSidebandCandidate * wRdiElectronSidebandCandidateFilter)
