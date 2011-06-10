import FWCore.ParameterSet.Config as cms

hNuGenFilter = cms.EDFilter (
    "HeavyNuGenFilter",
    keepIds = cms.vint32(1,)
    )
