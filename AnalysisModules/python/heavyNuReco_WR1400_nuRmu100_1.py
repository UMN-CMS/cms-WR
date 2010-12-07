import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1400_nuRmu100/HeavyNuGenHLT_WR1400_nuRmu100_1-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1400_nuRmu100/HeavyNuGenHLT_WR1400_nuRmu100_2-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1400_nuRmu100/HeavyNuGenHLT_WR1400_nuRmu100_3-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1400_nuRmu100/HeavyNuGenHLT_WR1400_nuRmu100_4-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1400_nuRmu100/HeavyNuGenHLT_WR1400_nuRmu100_5-reco-pool.root'
])
