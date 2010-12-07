import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu1200/HeavyNuGenHLT_WR1500_nuRmu1200_1-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu1200/HeavyNuGenHLT_WR1500_nuRmu1200_2-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu1200/HeavyNuGenHLT_WR1500_nuRmu1200_3-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu1200/HeavyNuGenHLT_WR1500_nuRmu1200_4-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu1200/HeavyNuGenHLT_WR1500_nuRmu1200_5-reco-pool.root'
])
