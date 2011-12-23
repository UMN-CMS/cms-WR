import FWCore.ParameterSet.Config as cms

hNuTop = cms.EDFilter(
    "HeavyNuTop",

    isPFJets = cms.bool( True ), 
    trigMatchPset = cms.PSet(
        trigEventTag = cms.InputTag( "" ),
        muonTriggers = cms.vstring( '' ),
        trigEra      = cms.int32( 20111 ),
        firstRun     = cms.vint32( 0 ),
        lastRun      = cms.vint32( 999999 ),
        muonMatch    = cms.string( '' ),
        triggerPt    = cms.double( 40. ),
        randomSeed   = cms.int32( 0 )  # for MC
    ),
    muIDPset = cms.PSet(
        eraForId     = cms.int32( 20111 )
    ),
    DoLog          = cms.bool( False ),
    muonTag        = cms.InputTag( 'selectedPatMuons' ),
    jetTag         = cms.InputTag( 'selectedPatJets' ),
    metTag         = cms.InputTag( 'patMETs' ),
    electronTag    = cms.InputTag( 'selectedPatElectrons' ),
    BtagName       = cms.string('trackCountingHighEffBJetTags'),
    minBtagDiscr   = cms.double(0.669), # yields 0.1% fake rate, see SWGuideBTagPerformance twiki
    minLep1pt      = cms.double(60.),
    minLep2pt      = cms.double(30.),
    minJetPt       = cms.double(40),
    maxMuAbsEta    = cms.double(2.4),
    maxElecAbsEta  = cms.double(2.5),
    maxJetAbsEta   = cms.double(2.5),
    minMuMuMass    = cms.double(200),
    min4objMass    = cms.double(520),
    minLeptonJetdR = cms.double(0.5),
    
    muonTrackRelIsoLimit  = cms.double(0.1), # 10.0),
    maxVertexZsepCM       = cms.double(0.03),
    maxVertexJetVZsepCM   = cms.double(0.1),

    pileupEra         = cms.int32(20111),

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
