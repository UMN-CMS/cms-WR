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
    fileNames=cms.untracked.vstring('file:/local/cms/phedex/store/mc/Summer11/WRToNuLeptonToLLJJ_MW-1000_MNu-600_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S4_START42_V11-v1/0000/2A7409FF-06B6-E011-8ADC-00215E222370.root')
#    fileNames=cms.untracked.vstring('file:/home/ugrad/pastika/cms/gen/powheg-hvq/POWHEG_PYTHIA6_ttbar_lnublnub_7TeV_cff_py_GEN.root')
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("genanal.root"),
)

process.load("HeavyNu.AnalysisModules.heavynugenlevel_cfi")

process.load("RecoJets.Configuration.RecoGenJets_cff")
process.load("RecoJets.Configuration.GenJetParticles_cff")

process.prepSeq=cms.Sequence(process.genParticlesForJetsNoMuNoNu*process.ak5GenJetsNoMuNoNu)

import HeavyNu.AnalysisModules.heavynugenlevel_cfi;


