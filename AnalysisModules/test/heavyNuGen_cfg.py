import FWCore.ParameterSet.Config as cms

import os

needGenParts=True

process = cms.Process("AGEN");

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
# source

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
     duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames=cms.untracked.vstring('file:/local/cms/user/pastika/powheg-hvq/fact1.0_ren1.0/powheg_pythia/powheg_pythia_005.root')
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("genanal.root"),
)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("HeavyNu.AnalysisModules.heavynugenlevel_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")

process.load("RecoJets.Configuration.RecoGenJets_cff")
process.load("RecoJets.Configuration.GenJetParticles_cff")

if (needGenParts):
    process.p = cms.Path( process.genParticles*process.genParticlesForJetsNoMuNoNu*process.ak5GenJetsNoMuNoNu*process.hNuGen )
else:
    process.p = cms.Path( process.genParticlesForJetsNoMuNoNu*process.ak5GenJetsNoMuNoNu*process.hNuGen )
