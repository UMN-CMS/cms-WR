import FWCore.ParameterSet.Config as cms

# Selection of electrons from the slimmedElectrons collection (PAT)
""" \addtogroup electronSkim_Group electronSkim sequences
@{
"""

wRpreselJets = cms.EDFilter("PATJetRefSelector",
		src = cms.InputTag("slimmedJets"),
		cut = cms.string("pt>20"),
		)

wRpreselJetFilter = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("wRpreselJets"),
                                 minNumber = cms.uint32(2)
                             )
