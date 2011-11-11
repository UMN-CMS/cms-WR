import FWCore.ParameterSet.Config as cms

hNuGen = cms.EDAnalyzer (
    "HeavyNuGenLevel",
    minMu1pt     = cms.double(60.),
    minMu2pt     = cms.double(30.),
    minJetPt     = cms.double(40),
    maxMuAbsEta  = cms.double(2.4),
    maxJetAbsEta = cms.double(2.5),
    minMuMuMass  = cms.double(200),
    min4objMass  = cms.double(520),
    minMuonJetdR = cms.double(0.5)    
    )
