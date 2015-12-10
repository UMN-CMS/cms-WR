import FWCore.ParameterSet.Config as cms

from ExoAnalysis.cmsWR.skimMuon_cff import *
from ExoAnalysis.cmsWR.skimElectron_cff import *
from ExoAnalysis.cmsWR.skimJets_cff import *

# Preselections
""" \addtogroup microAODSeq_Group Sequences to produce microAOD
@{
"""

#mixed flavour candidates
wReleMuCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                  decay = cms.string("wRleadingElectron wRsubleadingMuon"),
                                  role = cms.string("leading subleading"),
                                  checkCharge = cms.bool(False),
                                  # the cut on the pt of the daughter is to respect the order of leading and subleading:
                                  # if both electrons have pt>60 GeV there will be two di-electron candidates with inversed order
                                  cut = cms.string("mass > 0 && daughter(0).pt>daughter(1).pt"),
)

wRmuEleCandidate = cms.EDProducer("CandViewShallowCloneCombiner",
                                  decay = cms.string("wRleadingMuon wRsubleadingElectron"),
                                  role = cms.string("leading subleading"),
                                  checkCharge = cms.bool(False),
                                  # the cut on the pt of the daughter is to respect the order of leading and subleading:
                                  # if both electrons have pt>60 GeV there will be two di-electron candidates with inversed order
                                  cut = cms.string("mass > 0 && daughter(0).pt>daughter(1).pt"),
)

wRdiLeptonCandidate = cms.EDProducer("CandViewMerger",
                                     src = cms.VInputTag("wRdiElectronCandidate", "wRdiMuonCandidate", "wReleMuCandidate", "wRmuEleCandidate")
                                 )


# filter events with at least 2 leptons
wRdiLeptonFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("wRdiLeptonCandidate"),
                                minNumber = cms.uint32(1)
                            )


wRdiLeptonSkimSequence = cms.Sequence(
    (wRleadingElectron + wRleadingMuon + wRsubleadingElectron + wRsubleadingMuon) *
    (wReleMuCandidate + wRmuEleCandidate + wRdiElectronCandidate + wRdiMuonCandidate) *
    wRdiLeptonCandidate * wRdiLeptonFilter
)


wRdijetSkimSequence = cms.Sequence( wRpreselJets * wRpreselJetFilter )

