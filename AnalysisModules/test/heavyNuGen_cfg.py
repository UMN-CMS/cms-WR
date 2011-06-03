import FWCore.ParameterSet.Config as cms

import os

process = cms.Process("AGEN");

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
# source

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring('file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1000_nuRmu500/HeavyNuGenHLT_WR1000_nuRmu500_1-reco-pool.root')
#/hdfs/cms/skim/mu/39X/Dec22ReReco/Run2010B/Mu_Run2010B_Dec22ReReco_v1_AOD_068-muSkim-pool.root')
#This file is AT FNAL:
#    fileNames=cms.untracked.vstring('/store/mc/Fall10/Z1Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola/GEN-SIM-RECO/START38_V12-v2/0018/027050CA-D50B-E011-91DD-00261894391C.root')
#    fileNames=cms.untracked.vstring( "file:/hdfs/cms/user/heavynu/HeavyNuRecoFromHLT/38X/WR1000_nuRmu100/HeavyNuGenHLT_WR1000_nuRmu100_1-reco-pool.root" )
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("genanal.root"),
)

process.load("HeavyNu.AnalysisModules.heavynugenlevel_cfi")

process.load("RecoJets.Configuration.RecoGenJets_cff")
process.load("RecoJets.Configuration.GenJetParticles_cff")

process.p = cms.Path( process.genParticlesForJetsNoMuNoNu*process.ak5GenJetsNoMuNoNu*process.hNuGen )
