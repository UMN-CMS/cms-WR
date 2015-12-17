import FWCore.ParameterSet.Config as cms

from ExoAnalysis.cmsWR.skimMuon_cff import *
from ExoAnalysis.cmsWR.skimElectron_cff import *
from ExoAnalysis.cmsWR.skimJets_cff import *
from ExoAnalysis.cmsWR.skimEMu_cff import *

# Preselections
""" \addtogroup microAODSeq_Group Sequences to produce microAOD
@{
"""


wRdiLeptonCandidatePresel = cms.EDProducer("CandViewMerger",
                                     src = cms.VInputTag("wRdiElectronCandidatePresel", "wRdiMuonCandidatePresel", "wReleMuCandidatePresel", "wRmuEleCandidatePresel")
                                 )


# filter events with at least 2 leptons
wRdiLeptonFilterPresel = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("wRdiLeptonCandidatePresel"),
                                minNumber = cms.uint32(1)
                            )


wRdiLeptonSkimSequence = cms.Sequence(tunePMuons *
    (wRleadingElectronPresel + wRleadingMuonPresel + wRsubleadingElectronPresel + wRsubleadingMuonPresel) *
    (wReleMuCandidatePresel + wRmuEleCandidatePresel + wRdiElectronCandidatePresel + wRdiMuonCandidatePresel) *
    wRdiLeptonCandidatePresel * wRdiLeptonFilterPresel
)


wRdijetSkimSequence = cms.Sequence( wRpreselJets * wRpreselJetFilter )
