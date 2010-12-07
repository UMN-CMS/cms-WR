import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1600_nuRmu1500/HeavyNuGenHLT_WR1600_nuRmu1500_1-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1600_nuRmu1500/HeavyNuGenHLT_WR1600_nuRmu1500_2-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1600_nuRmu1500/HeavyNuGenHLT_WR1600_nuRmu1500_3-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1600_nuRmu1500/HeavyNuGenHLT_WR1600_nuRmu1500_4-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1600_nuRmu1500/HeavyNuGenHLT_WR1600_nuRmu1500_5-reco-pool.root'
])
