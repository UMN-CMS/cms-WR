import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1600_nuRmu700/HeavyNuGenHLT_WR1600_nuRmu700_1-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1600_nuRmu700/HeavyNuGenHLT_WR1600_nuRmu700_2-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1600_nuRmu700/HeavyNuGenHLT_WR1600_nuRmu700_3-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1600_nuRmu700/HeavyNuGenHLT_WR1600_nuRmu700_4-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1600_nuRmu700/HeavyNuGenHLT_WR1600_nuRmu700_5-reco-pool.root'
])
