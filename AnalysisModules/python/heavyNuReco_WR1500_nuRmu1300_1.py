import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu1300/HeavyNuGenHLT_WR1500_nuRmu1300_1-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu1300/HeavyNuGenHLT_WR1500_nuRmu1300_2-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu1300/HeavyNuGenHLT_WR1500_nuRmu1300_3-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu1300/HeavyNuGenHLT_WR1500_nuRmu1300_4-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu1300/HeavyNuGenHLT_WR1500_nuRmu1300_5-reco-pool.root'
])
