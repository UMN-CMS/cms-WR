import FWCore.ParameterSet.Config as cms

hNuTop = cms.EDFilter(
    "HeavyNuTop",
    trigMatchPset = cms.PSet(
        trigEventTag = cms.InputTag( "" ),
        muonMatch    = cms.string( '' ),
        randomSeed   = cms.int32( 0 ),  # for MC
        year         = cms.int32( 2011 ) # for MC
    ),
    DoLog        = cms.bool( False ),
    muonTag      = cms.InputTag( 'selectedPatMuons' ),
    jetTag       = cms.InputTag( 'selectedPatJets' ),
    metTag       = cms.InputTag( 'patMETs' ),
    electronTag  = cms.InputTag( 'selectedPatElectrons' ),
    BtagName     = cms.string('jetProbabilityBJetTags'),
    minBtagDiscr = cms.double(0.669), # yields 0.1% fake rate, see SWGuideBTagPerformance twiki
    minMu1pt     = cms.double(60.),
    minMu2pt     = cms.double(20.),
    minJetPt     = cms.double(40),
    maxMuAbsEta  = cms.double(2.4),
    maxEleAbsEta = cms.double(2.5),
    maxJetAbsEta = cms.double(2.5),
    minMuMuMass  = cms.double(200),
    min4objMass  = cms.double(520),
    minMuonJetdR = cms.double(0.5),
    
    muonTrackRelIsoLimit  = cms.double(0.1), # 10.0),
    maxVertexZsepCM       = cms.double(0.03),
    
    applyEleEScale    = cms.bool(False),
    EBscalefactor     = cms.double(1.0),
    EEscalefactor     = cms.double(1.0),
    applyEleIDweight  = cms.bool(False),
    EBidWgt           = cms.double(1.0),
    EEidWgt           = cms.double(1.0),
    applyMuIDEffcorr  = cms.bool(True),
    studyScaleFactor  = cms.bool(True),
    isSignal          = cms.bool(False),
    mNuRnormalization = cms.double(1000.0)
)
