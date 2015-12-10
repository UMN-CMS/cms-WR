import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi
wReejjHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*',
    ]
)


wRmumujjHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_Mu45_eta2p1_v*',
        'HLT_Mu50_v*',        
        ]
)

wRemujjHLTFilter =  HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v*',
        ]
)

tagAndProbeDoubleEleHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_Ele30WP60_Ele8_Mass55_v*',
        'HLT_Ele30WP60_SC4_Mass55_v*',
    ]
)

tagAndProbeDoubleMuHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_IsoMu22_v*',
        'HLT_IsoMu27_v*'
    ]
)


wRHLTFilterMC =  HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = wReejjHLTFilter.HLTPaths + wRmumujjHLTFilter.HLTPaths + wRemujjHLTFilter.HLTPaths
)

tagAndProbeHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = tagAndProbeDoubleMuHLTFilter.HLTPaths + tagAndProbeDoubleEleHLTFilter.HLTPaths
)
