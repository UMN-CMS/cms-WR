import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi


### \page triggers Triggers
## Set of triggers used in the analysis
## 
## 
## Here the list of triggers:
##  - \b HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*  
##    - electron channel
##  - \b HLT_Mu45_eta2p1_v* 
##    - actual muon channel and emjj sideband -> to be removed?
##  - \b HLT_Mu50_v*'                          
##    - [PROPOSED] new muon channel trigger
##  - \b HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v*  
##    - [PROPOSED] emujj sideband
##  - \b HLT_Ele30WP60_Ele8_Mass55_v*          
##    - electron tag&probe
##  - \b HLT_Ele30WP60_SC4_Mass55_v*           
##    -  electron tag&probe
##  - \b HLT_IsoMu22_v*                        
##    - [PROPOSED] muon tag&probe
##  - \b HLT_IsoMu27_v*                        
##    - [PROPOSED] muon tag&probe
## 

wReejjHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*', ### \ingroup hlt_Group electron channel trigger
    ]
)


wRmumujjHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_Mu45_eta2p1_v*', 
        'HLT_Mu50_v*',        ## \ingroup hlt_Group muon channel trigger
        ]
)

wRemujjHLTFilter =  HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v*', ## \ingroup hlt_Group flavour sideband trigger
        ]
)

tagAndProbeDoubleEleHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_Ele30WP60_Ele8_Mass55_v*',  ## \ingroup hlt_Group electron tagAndProbe trigger
        'HLT_Ele30WP60_SC4_Mass55_v*',
    ]
)

tagAndProbeDoubleMuHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_IsoMu22_v*', ## \ingroup hlt_Group muon tagAndProbe trigger
        'HLT_IsoMu27_v*' 
    ]
)


wRHLTFilter =  HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = wReejjHLTFilter.HLTPaths + wRmumujjHLTFilter.HLTPaths + wRemujjHLTFilter.HLTPaths
)

tagAndProbeHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = tagAndProbeDoubleMuHLTFilter.HLTPaths + tagAndProbeDoubleEleHLTFilter.HLTPaths
)


signalHltSequence = cms.Sequence(wRHLTFilter)
tagAndProbeHltSequence = cms.Sequence(tagAndProbeHLTFilter)

