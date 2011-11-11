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
    fileNames=cms.untracked.vstring('file:/home/ugrad/pastika/cms/gen/powheg-hvq/POWHEG_PYTHIA6_ttbar_lnublnub_7TeV_cff_py_GEN.root')
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("genanal.root"),
)

process.load("HeavyNu.AnalysisModules.heavynugenlevel_cfi")

process.load("RecoJets.Configuration.RecoGenJets_cff")
process.load("RecoJets.Configuration.GenJetParticles_cff")

process.p = cms.Path( process.genParticlesForJetsNoMuNoNu*process.ak5GenJetsNoMuNoNu*process.hNuGen )
