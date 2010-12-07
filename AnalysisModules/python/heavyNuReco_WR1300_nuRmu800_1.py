import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu800/HeavyNuGenHLT_WR1300_nuRmu800_1-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu800/HeavyNuGenHLT_WR1300_nuRmu800_2-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu800/HeavyNuGenHLT_WR1300_nuRmu800_3-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu800/HeavyNuGenHLT_WR1300_nuRmu800_4-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu800/HeavyNuGenHLT_WR1300_nuRmu800_5-reco-pool.root'
])
