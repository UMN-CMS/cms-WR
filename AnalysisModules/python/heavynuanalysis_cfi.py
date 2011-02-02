import FWCore.ParameterSet.Config as cms

hNu = cms.EDFilter(
    "HeavyNu",
    trigMatchPset = cms.PSet(
    trigEventTag = cms.InputTag( "" ),
    muonMatch    = cms.string( '' ),
    randomSeed   = cms.int32( 0 )  # for MC
    ),
    DoLog        = cms.bool( True ),
    muonTag      = cms.InputTag( 'selectedPatMuons' ),
    jetTag       = cms.InputTag( 'selectedPatJets' ),
    electronTag  = cms.InputTag( 'selectedPatElectrons' ),
    BtagName     = cms.string('jetProbabilityBJetTags'),
    minMu1pt     = cms.double(60.),
    minMu2pt     = cms.double(20.),
    minJetPt     = cms.double(40),
    maxMuAbsEta  = cms.double(2.4),
    maxJetAbsEta = cms.double(2.5),
    minMuMuMass  = cms.double(200),
    min4objMass  = cms.double(520),
    minMuonJetdR = cms.double(0.3), # (0.5)
    
    muonTrackIsoLimitGeV  = cms.double(10.0),
    maxVertexZsepCM       = cms.double(0.03),
    
    ZmassWinMinGeV= cms.double(84.),
    ZmassWinMaxGeV= cms.double(98.),

    applyJECUsign    = cms.int32(0),
    applyTrigEffsign = cms.int32(0),
    applyMESfactor   = cms.double(1.0),

    isSignal     = cms.bool(False),
    mNuRnormalization = cms.double(1000.0)
    )
