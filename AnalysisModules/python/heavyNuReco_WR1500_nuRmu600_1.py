import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu600/HeavyNuGenHLT_WR1500_nuRmu600_1-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu600/HeavyNuGenHLT_WR1500_nuRmu600_2-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu600/HeavyNuGenHLT_WR1500_nuRmu600_3-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu600/HeavyNuGenHLT_WR1500_nuRmu600_4-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1500_nuRmu600/HeavyNuGenHLT_WR1500_nuRmu600_5-reco-pool.root'
])
