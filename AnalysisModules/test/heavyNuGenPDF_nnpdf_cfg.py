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


process.pdf0 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf0.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf0.pdfReweightBaseId = 0
process.pdf0.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf0.pdfReweightTargetId = 0
process.pdf0.doPDFReweight = True

process.pdf1 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf1.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf1.pdfReweightBaseId = 0
process.pdf1.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf1.pdfReweightTargetId = 1
process.pdf1.doPDFReweight = True

process.pdf2 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf2.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf2.pdfReweightBaseId = 0
process.pdf2.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf2.pdfReweightTargetId = 2
process.pdf2.doPDFReweight = True

process.pdf3 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf3.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf3.pdfReweightBaseId = 0
process.pdf3.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf3.pdfReweightTargetId = 3
process.pdf3.doPDFReweight = True

process.pdf4 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf4.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf4.pdfReweightBaseId = 0
process.pdf4.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf4.pdfReweightTargetId = 4
process.pdf4.doPDFReweight = True

process.pdf5 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf5.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf5.pdfReweightBaseId = 0
process.pdf5.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf5.pdfReweightTargetId = 5
process.pdf5.doPDFReweight = True

process.pdf6 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf6.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf6.pdfReweightBaseId = 0
process.pdf6.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf6.pdfReweightTargetId = 6
process.pdf6.doPDFReweight = True

process.pdf7 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf7.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf7.pdfReweightBaseId = 0
process.pdf7.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf7.pdfReweightTargetId = 7
process.pdf7.doPDFReweight = True

process.pdf8 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf8.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf8.pdfReweightBaseId = 0
process.pdf8.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf8.pdfReweightTargetId = 8
process.pdf8.doPDFReweight = True

process.pdf9 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf9.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf9.pdfReweightBaseId = 0
process.pdf9.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf9.pdfReweightTargetId = 9
process.pdf9.doPDFReweight = True

process.pdf10 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf10.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf10.pdfReweightBaseId = 0
process.pdf10.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf10.pdfReweightTargetId = 10
process.pdf10.doPDFReweight = True

process.pdf11 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf11.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf11.pdfReweightBaseId = 0
process.pdf11.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf11.pdfReweightTargetId = 11
process.pdf11.doPDFReweight = True

process.pdf12 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf12.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf12.pdfReweightBaseId = 0
process.pdf12.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf12.pdfReweightTargetId = 12
process.pdf12.doPDFReweight = True

process.pdf13 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf13.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf13.pdfReweightBaseId = 0
process.pdf13.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf13.pdfReweightTargetId = 13
process.pdf13.doPDFReweight = True

process.pdf14 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf14.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf14.pdfReweightBaseId = 0
process.pdf14.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf14.pdfReweightTargetId = 14
process.pdf14.doPDFReweight = True

process.pdf15 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf15.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf15.pdfReweightBaseId = 0
process.pdf15.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf15.pdfReweightTargetId = 15
process.pdf15.doPDFReweight = True

process.pdf16 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf16.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf16.pdfReweightBaseId = 0
process.pdf16.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf16.pdfReweightTargetId = 16
process.pdf16.doPDFReweight = True

process.pdf17 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf17.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf17.pdfReweightBaseId = 0
process.pdf17.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf17.pdfReweightTargetId = 17
process.pdf17.doPDFReweight = True

process.pdf18 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf18.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf18.pdfReweightBaseId = 0
process.pdf18.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf18.pdfReweightTargetId = 18
process.pdf18.doPDFReweight = True

process.pdf19 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf19.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf19.pdfReweightBaseId = 0
process.pdf19.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf19.pdfReweightTargetId = 19
process.pdf19.doPDFReweight = True

process.pdf20 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf20.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf20.pdfReweightBaseId = 0
process.pdf20.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf20.pdfReweightTargetId = 20
process.pdf20.doPDFReweight = True

process.pdf21 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf21.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf21.pdfReweightBaseId = 0
process.pdf21.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf21.pdfReweightTargetId = 21
process.pdf21.doPDFReweight = True

process.pdf22 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf22.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf22.pdfReweightBaseId = 0
process.pdf22.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf22.pdfReweightTargetId = 22
process.pdf22.doPDFReweight = True

process.pdf23 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf23.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf23.pdfReweightBaseId = 0
process.pdf23.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf23.pdfReweightTargetId = 23
process.pdf23.doPDFReweight = True

process.pdf24 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf24.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf24.pdfReweightBaseId = 0
process.pdf24.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf24.pdfReweightTargetId = 24
process.pdf24.doPDFReweight = True

process.pdf25 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf25.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf25.pdfReweightBaseId = 0
process.pdf25.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf25.pdfReweightTargetId = 25
process.pdf25.doPDFReweight = True

process.pdf26 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf26.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf26.pdfReweightBaseId = 0
process.pdf26.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf26.pdfReweightTargetId = 26
process.pdf26.doPDFReweight = True

process.pdf27 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf27.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf27.pdfReweightBaseId = 0
process.pdf27.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf27.pdfReweightTargetId = 27
process.pdf27.doPDFReweight = True

process.pdf28 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf28.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf28.pdfReweightBaseId = 0
process.pdf28.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf28.pdfReweightTargetId = 28
process.pdf28.doPDFReweight = True

process.pdf29 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf29.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf29.pdfReweightBaseId = 0
process.pdf29.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf29.pdfReweightTargetId = 29
process.pdf29.doPDFReweight = True

process.pdf30 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf30.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf30.pdfReweightBaseId = 0
process.pdf30.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf30.pdfReweightTargetId = 30
process.pdf30.doPDFReweight = True

process.pdf31 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf31.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf31.pdfReweightBaseId = 0
process.pdf31.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf31.pdfReweightTargetId = 31
process.pdf31.doPDFReweight = True

process.pdf32 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf32.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf32.pdfReweightBaseId = 0
process.pdf32.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf32.pdfReweightTargetId = 32
process.pdf32.doPDFReweight = True

process.pdf33 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf33.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf33.pdfReweightBaseId = 0
process.pdf33.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf33.pdfReweightTargetId = 33
process.pdf33.doPDFReweight = True

process.pdf34 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf34.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf34.pdfReweightBaseId = 0
process.pdf34.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf34.pdfReweightTargetId = 34
process.pdf34.doPDFReweight = True

process.pdf35 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf35.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf35.pdfReweightBaseId = 0
process.pdf35.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf35.pdfReweightTargetId = 35
process.pdf35.doPDFReweight = True

process.pdf36 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf36.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf36.pdfReweightBaseId = 0
process.pdf36.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf36.pdfReweightTargetId = 36
process.pdf36.doPDFReweight = True

process.pdf37 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf37.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf37.pdfReweightBaseId = 0
process.pdf37.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf37.pdfReweightTargetId = 37
process.pdf37.doPDFReweight = True

process.pdf38 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf38.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf38.pdfReweightBaseId = 0
process.pdf38.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf38.pdfReweightTargetId = 38
process.pdf38.doPDFReweight = True

process.pdf39 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf39.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf39.pdfReweightBaseId = 0
process.pdf39.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf39.pdfReweightTargetId = 39
process.pdf39.doPDFReweight = True

process.pdf40 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf40.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf40.pdfReweightBaseId = 0
process.pdf40.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf40.pdfReweightTargetId = 40
process.pdf40.doPDFReweight = True

process.pdf41 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf41.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf41.pdfReweightBaseId = 0
process.pdf41.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf41.pdfReweightTargetId = 41
process.pdf41.doPDFReweight = True

process.pdf42 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf42.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf42.pdfReweightBaseId = 0
process.pdf42.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf42.pdfReweightTargetId = 42
process.pdf42.doPDFReweight = True

process.pdf43 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf43.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf43.pdfReweightBaseId = 0
process.pdf43.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf43.pdfReweightTargetId = 43
process.pdf43.doPDFReweight = True

process.pdf44 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf44.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf44.pdfReweightBaseId = 0
process.pdf44.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf44.pdfReweightTargetId = 44
process.pdf44.doPDFReweight = True

process.pdf45 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf45.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf45.pdfReweightBaseId = 0
process.pdf45.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf45.pdfReweightTargetId = 45
process.pdf45.doPDFReweight = True

process.pdf46 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf46.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf46.pdfReweightBaseId = 0
process.pdf46.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf46.pdfReweightTargetId = 46
process.pdf46.doPDFReweight = True

process.pdf47 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf47.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf47.pdfReweightBaseId = 0
process.pdf47.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf47.pdfReweightTargetId = 47
process.pdf47.doPDFReweight = True

process.pdf48 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf48.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf48.pdfReweightBaseId = 0
process.pdf48.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf48.pdfReweightTargetId = 48
process.pdf48.doPDFReweight = True

process.pdf49 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf49.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf49.pdfReweightBaseId = 0
process.pdf49.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf49.pdfReweightTargetId = 49
process.pdf49.doPDFReweight = True

process.pdf50 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf50.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf50.pdfReweightBaseId = 0
process.pdf50.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf50.pdfReweightTargetId = 50
process.pdf50.doPDFReweight = True

process.pdf51 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf51.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf51.pdfReweightBaseId = 0
process.pdf51.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf51.pdfReweightTargetId = 51
process.pdf51.doPDFReweight = True

process.pdf52 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf52.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf52.pdfReweightBaseId = 0
process.pdf52.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf52.pdfReweightTargetId = 52
process.pdf52.doPDFReweight = True

process.pdf53 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf53.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf53.pdfReweightBaseId = 0
process.pdf53.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf53.pdfReweightTargetId = 53
process.pdf53.doPDFReweight = True

process.pdf54 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf54.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf54.pdfReweightBaseId = 0
process.pdf54.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf54.pdfReweightTargetId = 54
process.pdf54.doPDFReweight = True

process.pdf55 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf55.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf55.pdfReweightBaseId = 0
process.pdf55.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf55.pdfReweightTargetId = 55
process.pdf55.doPDFReweight = True

process.pdf56 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf56.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf56.pdfReweightBaseId = 0
process.pdf56.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf56.pdfReweightTargetId = 56
process.pdf56.doPDFReweight = True

process.pdf57 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf57.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf57.pdfReweightBaseId = 0
process.pdf57.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf57.pdfReweightTargetId = 57
process.pdf57.doPDFReweight = True

process.pdf58 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf58.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf58.pdfReweightBaseId = 0
process.pdf58.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf58.pdfReweightTargetId = 58
process.pdf58.doPDFReweight = True

process.pdf59 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf59.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf59.pdfReweightBaseId = 0
process.pdf59.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf59.pdfReweightTargetId = 59
process.pdf59.doPDFReweight = True

process.pdf60 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf60.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf60.pdfReweightBaseId = 0
process.pdf60.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf60.pdfReweightTargetId = 60
process.pdf60.doPDFReweight = True

process.pdf61 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf61.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf61.pdfReweightBaseId = 0
process.pdf61.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf61.pdfReweightTargetId = 61
process.pdf61.doPDFReweight = True

process.pdf62 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf62.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf62.pdfReweightBaseId = 0
process.pdf62.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf62.pdfReweightTargetId = 62
process.pdf62.doPDFReweight = True

process.pdf63 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf63.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf63.pdfReweightBaseId = 0
process.pdf63.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf63.pdfReweightTargetId = 63
process.pdf63.doPDFReweight = True

process.pdf64 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf64.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf64.pdfReweightBaseId = 0
process.pdf64.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf64.pdfReweightTargetId = 64
process.pdf64.doPDFReweight = True

process.pdf65 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf65.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf65.pdfReweightBaseId = 0
process.pdf65.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf65.pdfReweightTargetId = 65
process.pdf65.doPDFReweight = True

process.pdf66 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf66.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf66.pdfReweightBaseId = 0
process.pdf66.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf66.pdfReweightTargetId = 66
process.pdf66.doPDFReweight = True

process.pdf67 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf67.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf67.pdfReweightBaseId = 0
process.pdf67.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf67.pdfReweightTargetId = 67
process.pdf67.doPDFReweight = True

process.pdf68 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf68.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf68.pdfReweightBaseId = 0
process.pdf68.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf68.pdfReweightTargetId = 68
process.pdf68.doPDFReweight = True

process.pdf69 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf69.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf69.pdfReweightBaseId = 0
process.pdf69.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf69.pdfReweightTargetId = 69
process.pdf69.doPDFReweight = True

process.pdf70 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf70.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf70.pdfReweightBaseId = 0
process.pdf70.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf70.pdfReweightTargetId = 70
process.pdf70.doPDFReweight = True

process.pdf71 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf71.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf71.pdfReweightBaseId = 0
process.pdf71.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf71.pdfReweightTargetId = 71
process.pdf71.doPDFReweight = True

process.pdf72 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf72.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf72.pdfReweightBaseId = 0
process.pdf72.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf72.pdfReweightTargetId = 72
process.pdf72.doPDFReweight = True

process.pdf73 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf73.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf73.pdfReweightBaseId = 0
process.pdf73.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf73.pdfReweightTargetId = 73
process.pdf73.doPDFReweight = True

process.pdf74 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf74.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf74.pdfReweightBaseId = 0
process.pdf74.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf74.pdfReweightTargetId = 74
process.pdf74.doPDFReweight = True

process.pdf75 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf75.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf75.pdfReweightBaseId = 0
process.pdf75.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf75.pdfReweightTargetId = 75
process.pdf75.doPDFReweight = True

process.pdf76 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf76.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf76.pdfReweightBaseId = 0
process.pdf76.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf76.pdfReweightTargetId = 76
process.pdf76.doPDFReweight = True

process.pdf77 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf77.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf77.pdfReweightBaseId = 0
process.pdf77.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf77.pdfReweightTargetId = 77
process.pdf77.doPDFReweight = True

process.pdf78 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf78.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf78.pdfReweightBaseId = 0
process.pdf78.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf78.pdfReweightTargetId = 78
process.pdf78.doPDFReweight = True

process.pdf79 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf79.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf79.pdfReweightBaseId = 0
process.pdf79.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf79.pdfReweightTargetId = 79
process.pdf79.doPDFReweight = True

process.pdf80 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf80.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf80.pdfReweightBaseId = 0
process.pdf80.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf80.pdfReweightTargetId = 80
process.pdf80.doPDFReweight = True

process.pdf81 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf81.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf81.pdfReweightBaseId = 0
process.pdf81.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf81.pdfReweightTargetId = 81
process.pdf81.doPDFReweight = True

process.pdf82 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf82.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf82.pdfReweightBaseId = 0
process.pdf82.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf82.pdfReweightTargetId = 82
process.pdf82.doPDFReweight = True

process.pdf83 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf83.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf83.pdfReweightBaseId = 0
process.pdf83.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf83.pdfReweightTargetId = 83
process.pdf83.doPDFReweight = True

process.pdf84 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf84.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf84.pdfReweightBaseId = 0
process.pdf84.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf84.pdfReweightTargetId = 84
process.pdf84.doPDFReweight = True

process.pdf85 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf85.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf85.pdfReweightBaseId = 0
process.pdf85.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf85.pdfReweightTargetId = 85
process.pdf85.doPDFReweight = True

process.pdf86 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf86.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf86.pdfReweightBaseId = 0
process.pdf86.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf86.pdfReweightTargetId = 86
process.pdf86.doPDFReweight = True

process.pdf87 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf87.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf87.pdfReweightBaseId = 0
process.pdf87.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf87.pdfReweightTargetId = 87
process.pdf87.doPDFReweight = True

process.pdf88 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf88.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf88.pdfReweightBaseId = 0
process.pdf88.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf88.pdfReweightTargetId = 88
process.pdf88.doPDFReweight = True

process.pdf89 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf89.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf89.pdfReweightBaseId = 0
process.pdf89.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf89.pdfReweightTargetId = 89
process.pdf89.doPDFReweight = True

process.pdf90 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf90.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf90.pdfReweightBaseId = 0
process.pdf90.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf90.pdfReweightTargetId = 90
process.pdf90.doPDFReweight = True

process.pdf91 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf91.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf91.pdfReweightBaseId = 0
process.pdf91.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf91.pdfReweightTargetId = 91
process.pdf91.doPDFReweight = True

process.pdf92 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf92.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf92.pdfReweightBaseId = 0
process.pdf92.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf92.pdfReweightTargetId = 92
process.pdf92.doPDFReweight = True

process.pdf93 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf93.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf93.pdfReweightBaseId = 0
process.pdf93.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf93.pdfReweightTargetId = 93
process.pdf93.doPDFReweight = True

process.pdf94 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf94.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf94.pdfReweightBaseId = 0
process.pdf94.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf94.pdfReweightTargetId = 94
process.pdf94.doPDFReweight = True

process.pdf95 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf95.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf95.pdfReweightBaseId = 0
process.pdf95.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf95.pdfReweightTargetId = 95
process.pdf95.doPDFReweight = True

process.pdf96 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf96.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf96.pdfReweightBaseId = 0
process.pdf96.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf96.pdfReweightTargetId = 96
process.pdf96.doPDFReweight = True

process.pdf97 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf97.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf97.pdfReweightBaseId = 0
process.pdf97.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf97.pdfReweightTargetId = 97
process.pdf97.doPDFReweight = True

process.pdf98 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf98.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf98.pdfReweightBaseId = 0
process.pdf98.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf98.pdfReweightTargetId = 98
process.pdf98.doPDFReweight = True

process.pdf99 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf99.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf99.pdfReweightBaseId = 0
process.pdf99.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf99.pdfReweightTargetId = 99
process.pdf99.doPDFReweight = True

process.pdf100 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf100.pdfReweightBaseName = 'NNPDF20_100.LHgrid'
process.pdf100.pdfReweightBaseId = 0
process.pdf100.pdfReweightTargetName = 'NNPDF20_100.LHgrid'
process.pdf100.pdfReweightTargetId = 100
process.pdf100.doPDFReweight = True

process.p=cms.Path(process.prepSeq
     + process.pdf0
     + process.pdf1
     + process.pdf2
     + process.pdf3
     + process.pdf4
     + process.pdf5
     + process.pdf6
     + process.pdf7
     + process.pdf8
     + process.pdf9
     + process.pdf10
     + process.pdf11
     + process.pdf12
     + process.pdf13
     + process.pdf14
     + process.pdf15
     + process.pdf16
     + process.pdf17
     + process.pdf18
     + process.pdf19
     + process.pdf20
     + process.pdf21
     + process.pdf22
     + process.pdf23
     + process.pdf24
     + process.pdf25
     + process.pdf26
     + process.pdf27
     + process.pdf28
     + process.pdf29
     + process.pdf30
     + process.pdf31
     + process.pdf32
     + process.pdf33
     + process.pdf34
     + process.pdf35
     + process.pdf36
     + process.pdf37
     + process.pdf38
     + process.pdf39
     + process.pdf40
     + process.pdf41
     + process.pdf42
     + process.pdf43
     + process.pdf44
     + process.pdf45
     + process.pdf46
     + process.pdf47
     + process.pdf48
     + process.pdf49
     + process.pdf50
     + process.pdf51
     + process.pdf52
     + process.pdf53
     + process.pdf54
     + process.pdf55
     + process.pdf56
     + process.pdf57
     + process.pdf58
     + process.pdf59
     + process.pdf60
     + process.pdf61
     + process.pdf62
     + process.pdf63
     + process.pdf64
     + process.pdf65
     + process.pdf66
     + process.pdf67
     + process.pdf68
     + process.pdf69
     + process.pdf70
     + process.pdf71
     + process.pdf72
     + process.pdf73
     + process.pdf74
     + process.pdf75
     + process.pdf76
     + process.pdf77
     + process.pdf78
     + process.pdf79
     + process.pdf80
     + process.pdf81
     + process.pdf82
     + process.pdf83
     + process.pdf84
     + process.pdf85
     + process.pdf86
     + process.pdf87
     + process.pdf88
     + process.pdf89
     + process.pdf90
     + process.pdf91
     + process.pdf92
     + process.pdf93
     + process.pdf94
     + process.pdf95
     + process.pdf96
     + process.pdf97
     + process.pdf98
     + process.pdf99
     + process.pdf100
  )

