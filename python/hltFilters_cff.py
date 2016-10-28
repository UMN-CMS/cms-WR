import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt


### \page triggers Triggers
## Set of triggers used in the analysis
## 
## 
## Here the list of triggers:
##  - \b HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*  
##    - electron channel
##  - \b HLT_Mu50_v*'                          
##    - muon channel trigger
##  - \b HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v*  
##    - emujj sideband
##  - \b HLT_Ele30WP60_Ele8_Mass55_v*          
##    - electron tag&probe
##  - \b HLT_Ele30WP60_SC4_Mass55_v*           
##    -  electron tag&probe
##  - \b HLT_IsoMu22_v*                        
##    - muon tag&probe
##  - \b HLT_IsoMu27_v*                        
##    - muon tag&probe
## 

wReejjHLTFilter = hlt.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*', ### \ingroup hlt_Group electron channel trigger
    ]
)


wRmumujjHLTFilter = hlt.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_Mu50_v*',        ## \ingroup hlt_Group muon channel trigger
        ]
)

wRemujjHLTFilter =  hlt.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v*', ## \ingroup hlt_Group flavour sideband trigger
        ]
)

tagAndProbeDoubleEleHLTFilter = hlt.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_Ele*WP60_Ele8_Mass55_v*',  ## \ingroup hlt_Group electron tagAndProbe trigger
        'HLT_Ele*WP60_SC4_Mass55_v*',
    ]
)

tagAndProbeDoubleMuHLTFilter = hlt.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        'HLT_IsoMu22_v*', ## \ingroup hlt_Group muon tagAndProbe trigger
        'HLT_IsoMu27_v*' 
    ]
)


wRHLTFilter_MC =  hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT2"),
    throw = cms.bool(False),
    HLTPaths = wReejjHLTFilter.HLTPaths + wRmumujjHLTFilter.HLTPaths + wRemujjHLTFilter.HLTPaths
)

tagAndProbeHLTFilter_MC = hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT2"),
    throw = cms.bool(False),
    HLTPaths = tagAndProbeDoubleMuHLTFilter.HLTPaths + tagAndProbeDoubleEleHLTFilter.HLTPaths
)

wRHLTFilter_data =  hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    throw = cms.bool(False),
    HLTPaths = wReejjHLTFilter.HLTPaths + wRmumujjHLTFilter.HLTPaths + wRemujjHLTFilter.HLTPaths
)

tagAndProbeHLTFilter_data = hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    throw = cms.bool(False),
    HLTPaths = tagAndProbeDoubleMuHLTFilter.HLTPaths + tagAndProbeDoubleEleHLTFilter.HLTPaths
)

