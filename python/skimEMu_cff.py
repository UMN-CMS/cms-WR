import FWCore.ParameterSet.Config as cms
from ExoAnalysis.cmsWR.skimMuon_cff import *
from ExoAnalysis.cmsWR.skimElectron_cff import *

# Selection of muons from the slimmedMuons collection (PAT)
""" \addtogroup elemuonSkim_Group elemuonSkim sequences
@{
"""

### create ele-muon pair in signal region
wREleMuonCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                       decay = cms.string("wRsubleadingMuon wRsubleadingElectron"),
                                       role = cms.string("leading subleading"),
                                       checkCharge = cms.bool(False),
                                       # the cut on the pt of the daughter is to respect the order of leading and subleading:
                                           # if both muons have pt>60 GeV there will be two di-muon candidates with inversed order
                                       cut = cms.string("mass > 200 && daughter(0).pt>daughter(1).pt"),
                                       )

### filter: at least one ele-muon candidate in signal region
wREleMuonCandidateFilter = cms.EDFilter("CandViewCountFilter",
                                           src = cms.InputTag("wREleMuonCandidate"),
                                           minNumber = cms.uint32(1)
                                           )

### create ele-muon pair in sideband region
wREleMuonSidebandCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                       decay = cms.string("wRsubleadingMuon wRsubleadingElectron"),
                                       checkCharge = cms.bool(False),
                                       cut = cms.string("0< mass < 200 && daughter(0).pt>daughter(1).pt")
                                       )

### filter: at least one ele-muon candidate in sideband region
wREleMuonSidebandCandidateFilter = cms.EDFilter("CandViewCountFilter",
                                                   src = cms.InputTag("wREleMuonSidebandCandidate"),
                                                   minNumber = cms.uint32(1)
                                                   )

### di-muon selection sequence
wRmuonSelectionSeq = cms.Sequence(wRtunePMuons + wRleadingMuon + wRsubleadingMuon + wRlooseJet)
### di-electron selection sequence
wRelectronSelectionSeq = cms.Sequence(wRleadingElectron + wRsubleadingElectron)
### di-muon selection in signal region sequence
wREleMuonSignalSeq = cms.Sequence((wRmuonSelectionSeq + wRelectronSelectionSeq) * wREleMuonCandidate * wREleMuonCandidateFilter)
### di-muon selection in sideband region sequence
wREleMuonSidebandSeq = cms.Sequence((wRmuonSelectionSeq + wRelectronSelectionSeq) * wREleMuonSidebandCandidate * wREleMuonSidebandCandidateFilter)


### @}
