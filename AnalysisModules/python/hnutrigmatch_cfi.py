import FWCore.ParameterSet.Config as cms

# PAT trigger
muonTriggerMatchHLTMuons = cms.EDProducer (
    "PATTriggerMatcherDRDPtLessByR",
    src            = cms.InputTag( 'selectedPatMuons' ),
    matched        = cms.InputTag( 'patTrigger' ),
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( 'TriggerMuon' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring(),
    collectionTags = cms.vstring( '*' ),
    maxDPtRel      = cms.double( 1.0 ),
    maxDeltaR      = cms.double( 0.2 ),
    maxDeltaEta    = cms.double( 0.2 ), # no effect here
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True )
)
