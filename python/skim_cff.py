import FWCore.ParameterSet.Config as cms

from ExoAnalysis.cmsWR.skimMuon_cff import *
from ExoAnalysis.cmsWR.skimElectron_cff import *
from ExoAnalysis.cmsWR.skimJets_cff import *
from ExoAnalysis.cmsWR.skimEMu_cff import *

# Preselections
""" \addtogroup microAODSeq_Group Sequences to produce microAOD
@{
"""


wRdiLeptonCandidate = cms.EDProducer("CandViewMerger",
                                     src = cms.VInputTag("wRdiElectronCandidate", "wRdiMuonCandidate", "wReleMuCandidate", "wRmuEleCandidate")
                                 )


# filter events with at least 2 leptons
wRdiLeptonFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("wRdiLeptonCandidate"),
                                minNumber = cms.uint32(1)
                            )


wRdiLeptonSkimSequence = cms.Sequence(tunePMuons *
    (wRleadingElectron + wRleadingMuon + wRsubleadingElectron + wRsubleadingMuon) *
    (wReleMuCandidate + wRmuEleCandidate + wRdiElectronCandidate + wRdiMuonCandidate) *
    wRdiLeptonCandidate * wRdiLeptonFilter
)


wRdijetSkimSequence = cms.Sequence( wRpreselJets * wRpreselJetFilter )
