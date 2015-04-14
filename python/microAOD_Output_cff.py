import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff import MINIAODSIMEventContent

MICROAODSIMEventContent= cms.PSet(
    # put this if you have a filter
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('')
        ),
    outputCommands = MINIAODSIMEventContent.outputCommands+
    [
        'drop HcalNoiseSummary_*_*_*',
        'drop patJets_slimmedJetsAK8PFCHSSoftDropPacked_*_*',
        'drop patJets_slimmedJetsCMSTopTagCHSPacked_*_*',
        'drop patJets_slimmedJetsAK8_*_*',
        'drop patTaus_slimmedTaus_*_*',
        'drop patPackedCandidates_packedPFCandidates_*_*',
        'drop patPackedCandidates_lostTracks_*_*',
        'drop recoGenJets_slimmedGenJetsAK8_*_*',
        'drop recoCATopJetTagInfos_caTopTagInfosPAT_*_*',
        'drop l1extraL1HFRingss_l1extraParticles_*_*',
        ]
    )

