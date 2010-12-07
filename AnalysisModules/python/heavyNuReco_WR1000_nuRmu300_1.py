import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu300/HeavyNuGenHLT_WR1000_nuRmu300_1-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu300/HeavyNuGenHLT_WR1000_nuRmu300_2-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu300/HeavyNuGenHLT_WR1000_nuRmu300_3-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu300/HeavyNuGenHLT_WR1000_nuRmu300_4-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu300/HeavyNuGenHLT_WR1000_nuRmu300_5-reco-pool.root'
])
