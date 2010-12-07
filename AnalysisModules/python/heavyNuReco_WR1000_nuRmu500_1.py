import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu500/HeavyNuGenHLT_WR1000_nuRmu500_1-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu500/HeavyNuGenHLT_WR1000_nuRmu500_2-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu500/HeavyNuGenHLT_WR1000_nuRmu500_3-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu500/HeavyNuGenHLT_WR1000_nuRmu500_4-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu500/HeavyNuGenHLT_WR1000_nuRmu500_5-reco-pool.root'
])
