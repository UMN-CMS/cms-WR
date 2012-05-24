import FWCore.ParameterSet.Config as cms

HeavyNuEleTriggerEff = cms.EDAnalyzer(
    'HeavyNuEleTriggerEff',
    
    trigEventTag     = cms.InputTag("patTriggerEvent"),
    doDebugMessages  = cms.bool(False),
    #doDebugMessages  = cms.bool(True),
    
    # to implement HEEP ID
    # rhoForHEEPId      = cms.InputTag("kt6PFJetsForIsolation","rho","triggerEfficienciesEle"),   # this worked in the standalone version
    #rhoForHEEPId      = cms.InputTag( 'kt6PFJetsForIsolation','rho' ),                            # this should be ok for the usage integrated within heavyNuAnalysis_cfg.py 
    #rhoForHEEPId      = cms.InputTag( 'selectedPatJetsPFlow','rho','PAT' ),                           # this should be ok for the usage integrated within heavyNuAnalysis_cfg.py 
    rhoForHEEPId      = cms.InputTag( "kt6PFJetsForIsolation", "rho", "PAT" ),                           # this should be ok for the usage integrated within heavyNuAnalysis_cfg.py 
    
    electronTag  = cms.InputTag("patElectrons"),
    maxAbsEtaOfflEle = cms.double( 2.5 ),
    minPtOfflEle     = cms.double(  40 ),
    heepVersion      = cms.int32(   40 ),
    ebScale          = cms.double(   1. ),
    eeScale          = cms.double(   1. ),

    #jetTag           = cms.InputTag("patJets","pfCandidates",""),
    #jetTag           = cms.InputTag("cleanPatJets","",""), # this worked in the standalone version
    #jetTag           = cms.InputTag("kt6PFJetsForIsolation"),
    #jetTag           = cms.InputTag("kt6PFJetsForIsolation","","PAT"),
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

    #electronFilters = cms.vstring('hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter')
    # the last filter of 'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_vx'
    #electronSeedingFilters = cms.vstring('HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17Mass50Sequence'),
    #electronSeedingFilters = cms.vstring('HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17Mass50Sequence'),


    electronSeedingFilters = cms.vstring('hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter'),
    # documentation of this choice partly addressed here: 
    # https://hypernews.cern.ch/HyperNews/CMS/get/egamma-hlt/168/1.html

    
    electronTargetFilters  = cms.vstring('hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter'),
  
    
    plotFolderName  = cms.string('twoJets'),
    
    )
