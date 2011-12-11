import FWCore.ParameterSet.Config as cms

hNuQCD = cms.EDFilter(
    "MuJetBackground",
    DoLog        = cms.bool( False ),
    trigMatchPset = cms.PSet(
        trigEventTag = cms.InputTag( "" ),
        muonTriggers = cms.vstring( '' ),
        triggerPt    = cms.double( 40. ),
        trigEra      = cms.int32( 20110 ),
        firstRun     = cms.vint32( 0 ),
        lastRun      = cms.vint32( 999999 ),
        muonMatch    = cms.string( '' ),
        randomSeed   = cms.int32( 0 )  # for MC
    ),
    muIDPset = cms.PSet(
        eraForId     = cms.int32( 2011 )
    ),
    pileupEra    = cms.int32(20119),
    isPFJets     = cms.bool(True),
    muonTag      = cms.InputTag( 'selectedPatMuons' ),
    jetTag       = cms.InputTag( 'selectedPatJets' ),
    jptTag       = cms.InputTag( 'JetPlusTrackZSPCorJetAntiKt5' ),
    minMu1pt     = cms.double(60.),
    minMu2pt     = cms.double(30.),
    minJetPt     = cms.double(40),
    maxMuAbsEta  = cms.double(2.4),
    maxJetAbsEta = cms.double(2.5),
    minMuMuMass  = cms.double(200),
    minMuonJetdR = cms.double(0.5),
    
    muonTrackRelIsoLimit  = cms.double(0.1), 
    dimuonMaxVertexZsepCM = cms.double(0.03),
    maxJetVZsepCM         = cms.double(0.1),
    
    #--- Specific variables added for INR QCD cross-check ---#
    electronTag = cms.InputTag( 'selectedPatElectrons' ),
    # photonTag   = cms.InputTag( 'selectedPatElectrons' ),
    # hybridSCs   = cms.InputTag( 'correctedHybridSuperClusters' ),
    # multi5x5SCs = cms.InputTag( 'correctedMulti5x5SuperClustersWithPreshower' ),
    minimumSuperClusterEt = cms.double(10.0), 

    getSurvivalRate = cms.bool(False),    
    doClosureTest   = cms.bool(False),    
    doQuadJetTest   = cms.bool(False),    

    minimumMuJetdPhi          = cms.double(2.8274334),
    minimumJetPtForDijets     = cms.double(20.),
    minimumDeltaRforExtraJets = cms.double(0.7),

    reweightPtLow  = cms.vdouble( 30.0 ),
    reweightPtHigh = cms.vdouble( 10000.0 ),
    reweightTight  = cms.vdouble( 1.0 ),

    # Take PF MET 
    METvariety = cms.int32(2) 
)
