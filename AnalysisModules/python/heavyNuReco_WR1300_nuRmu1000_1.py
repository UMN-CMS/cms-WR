import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu1000/HeavyNuGenHLT_WR1300_nuRmu1000_1-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu1000/HeavyNuGenHLT_WR1300_nuRmu1000_2-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu1000/HeavyNuGenHLT_WR1300_nuRmu1000_3-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu1000/HeavyNuGenHLT_WR1300_nuRmu1000_4-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu1000/HeavyNuGenHLT_WR1300_nuRmu1000_5-reco-pool.root'
])
