import FWCore.ParameterSet.Config as cms

HeavyNuEleTriggerEff = cms.EDAnalyzer(
    'HeavyNuEleTriggerEff',


    runMin      = cms.int32(   0 ),
    #runMax      = cms.int32(   191717),
    #runMin      = cms.int32(   191718),
    runMax      = cms.int32(   999999999 ),

    # this is the HLT object
    trigEventHLTTag     = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    # this is the pat object
    trigEventPatTag     = cms.InputTag("patTriggerEvent"),

    doDebugMessages  = cms.bool(False),
    #doDebugMessages  = cms.bool(True),
    
    # to implement HEEP ID
    rhoForHEEPId      = cms.InputTag( "kt6PFJetsForIsolation", "rho", "PAT" ),
    
    electronTag  = cms.InputTag("patElectrons"),
    maxAbsEtaOfflEle = cms.double( 2.5 ),
    minPtOfflEle     = cms.double(  40 ),

    # 40 is the HEEP version customary for signal selection of the 2012 analysis 
    heepVersion      = cms.int32(   40 ),
    # 0 is meant to avoid any HEEP selection, just leave gsf electrons as they come from pat 
    #heepVersion      = cms.int32(    0 ),
    
    ebScale          = cms.double(   1. ),
    eeScale          = cms.double(   1. ),

    jetTag           = cms.InputTag("selectedPatJetsPFlow", "","PAT"),
    numOfflJets      = cms.int32(    2 ),
    minPtOfflJets    = cms.double(   40),
    
    
    seedHLTForL1Efficiency =  cms.vstring(
    'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v3',
    'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v4',
    'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v5'
    ),

    targetL1Algo          = cms.string('L1_SingleEG22'),

    nonIsolatedEmSource = cms.InputTag("l1extraParticles","NonIsolated"),
    isolatedEmSource    = cms.InputTag("l1extraParticles","Isolated"),
    
    seedHLTForHLTEfficiency =  cms.vstring(
    'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v3',
    'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v4',
    'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v5'
    ),
    minPtOfObjects = cms.double( 33. ),

    targetHLTPaths          = cms.vstring(
    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v2',
    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3',
    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v4', # this seems to be the matchging one
    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v5',
    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v6',
    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7',
    ),


    # this filter is the last filter in the path:
    #                                            HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v3
    electronSeedingPathEndFilter = cms.vstring('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter'),
    # documentation of this choice partly addressed here: 
    # https://hypernews.cern.ch/HyperNews/CMS/get/egamma-hlt/168/1.html
    #electronSeedingFilters = cms.vstring('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter'),

    # filter which is fired by the 'good electron' of the tag-and-probe seeding HLT path
    electronSeedingPathTagFilter = cms.vstring('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter'),

    
    electronTargetFilters  = cms.vstring('hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter'),
  
    
    plotFolderName  = cms.string('twoJets'),
    
    )
