import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/jmmans/PYTHIA6_Tauola_TTbar_mumu_TuneZ2_7TeV_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RECO_385/PYTHIA6_Tauola_TTbar_mumu_TuneZ2_7TeV_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RECO_118.root'
])
