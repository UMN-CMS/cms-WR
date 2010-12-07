import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1200_nuRmu900/HeavyNuGenHLT_WR1200_nuRmu900_1-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1200_nuRmu900/HeavyNuGenHLT_WR1200_nuRmu900_2-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1200_nuRmu900/HeavyNuGenHLT_WR1200_nuRmu900_3-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1200_nuRmu900/HeavyNuGenHLT_WR1200_nuRmu900_4-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1200_nuRmu900/HeavyNuGenHLT_WR1200_nuRmu900_5-reco-pool.root'
])
