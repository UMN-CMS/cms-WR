import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1400_nuRmu200/HeavyNuGenHLT_WR1400_nuRmu200_1-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1400_nuRmu200/HeavyNuGenHLT_WR1400_nuRmu200_2-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1400_nuRmu200/HeavyNuGenHLT_WR1400_nuRmu200_3-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1400_nuRmu200/HeavyNuGenHLT_WR1400_nuRmu200_4-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1400_nuRmu200/HeavyNuGenHLT_WR1400_nuRmu200_5-reco-pool.root'
])
