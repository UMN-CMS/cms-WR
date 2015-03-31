import FWCore.ParameterSet.Config as cms

# Selection of muons from the slimmedMuons collection (PAT)
""" \addtogroup muonSkim_Group muonSkim sequences
@{
"""

## create new muons with the tuneP track \ingroup muonSkim_Group
wRtunePMuons = cms.EDProducer("TunePMuonProducer",
                             src = cms.InputTag("slimmedMuons")
                             )


### select leading muon
wRleadingMuon = cms.EDFilter("CandViewSelector",
                                 src = cms.InputTag("wRtunePMuons"),
                                 cut = cms.string("pt>60"),
                                 )

### select subleading muon
wRsubleadingMuon = cms.EDFilter("CandViewSelector",
                                 src = cms.InputTag("wRtunePMuons"),
                                 cut = cms.string("pt>40"),
                                 )
### create di-muon pair in signal region
wRdiMuonCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                       decay = cms.string("wRleadingMuon wRsubleadingMuon"),
                                       role = cms.string("leading subleading"),
                                       checkCharge = cms.bool(False),
                                       # the cut on the pt of the daughter is to respect the order of leading and subleading:
                                           # if both muons have pt>60 GeV there will be two di-muon candidates with inversed order
                                       cut = cms.string("mass > 200 && daughter(0).pt>daughter(1).pt"),
                                       )

### filter: at least one di-muon candidate in signal region
wRdiMuonCandidateFilter = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("wRdiMuonCandidate"),
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


### di-muon selection sequence
wRmuonSelectionSeq = cms.Sequence(wRleadingMuon + wRsubleadingMuon)
### di-muon selection in signal region sequence
wRdiMuonSignalSeq = cms.Sequence(wRmuonSelectionSeq * wRdiMuonCandidate * wRdiMuonCandidateFilter)
### di-muon selection in sideband region sequence
wRdiMuonSidebandSeq = cms.Sequence(wRmuonSelectionSeq * wRdiMuonSidebandCandidate * wRdiMuonSidebandCandidateFilter)


### @}
