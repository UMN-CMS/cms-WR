import FWCore.ParameterSet.Config as cms

hNu = cms.EDFilter(
    "HeavyNu",
    trigMatchPset = cms.PSet(
        trigEventTag = cms.InputTag( "" ),
        muonTriggers = cms.vstring( '' ),
        firstRun     = cms.vint32( 0 ),
        lastRun      = cms.vint32( 999999 ),
        muonMatch    = cms.string( '' ),
        triggerPt    = cms.double( 40. ),
        randomSeed   = cms.int32( 0 ),  # for MC
        year         = cms.int32( 2011 ) # for MC
    ),
    muIDPset = cms.PSet(
        eraForId     = cms.int32( 2011 )
    ),
    DoLog        = cms.bool( False ),
    muonTag      = cms.InputTag( 'selectedPatMuons' ),
    jetTag       = cms.InputTag( 'selectedPatJets' ),
    metTag       = cms.InputTag( 'patMETs' ),
    electronTag  = cms.InputTag( 'selectedPatElectrons' ),
    trackTag     = cms.InputTag( 'patTracksPt10' ),
    BtagName     = cms.string('jetProbabilityBJetTags'),
    minBtagDiscr = cms.double(0.669), # yields 0.1% fake rate, see SWGuideBTagPerformance twiki
    minMu1pt     = cms.double(60.),
    minMu2pt     = cms.double(20.),
    minJetPt     = cms.double(40),
    maxMuAbsEta  = cms.double(2.4),
    maxJetAbsEta = cms.double(2.5),
    minMuMuMass  = cms.double(200),
    min4objMass  = cms.double(520),
    minMuonJetdR = cms.double(0.5),
    
    muonTrackRelIsoLimit  = cms.double(0.1), # 10.0),
    maxVertexZsepCM       = cms.double(0.03),
    maxJetVZsepCM         = cms.double(0.1),
    
    ZmassWinMinGeV= cms.double(86.),
    ZmassWinMaxGeV= cms.double(96.),

    pileupEra         = cms.int32(20110),
    DisableTriggerCorrection = cms.bool(False),

    jecEra            = cms.int32(0),
    applyJECUsign     = cms.int32(0),
    applyTrigEffsign  = cms.int32(0),
    studyMuSelectEff  = cms.bool(True),
    applyMESfactor    = cms.double(1.0),
    applyMuIDEffcorr  = cms.bool(False),
    applyMuIDEffsign  = cms.int32(0),

    studyScaleFactor  = cms.bool(True),
    highestPtTriggerOnly = cms.bool(False),
    isSignal     = cms.bool(False),
    mNuRnormalization = cms.double(1000.0),

    isPFJets = cms.bool(False)
    )
