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
process.pdf0.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf0.pdfReweightBaseId = 0
process.pdf0.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf0.pdfReweightTargetId = 0
process.pdf0.doPDFReweight = True

process.pdf1 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf1.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf1.pdfReweightBaseId = 0
process.pdf1.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf1.pdfReweightTargetId = 1
process.pdf1.doPDFReweight = True

process.pdf2 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf2.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf2.pdfReweightBaseId = 0
process.pdf2.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf2.pdfReweightTargetId = 2
process.pdf2.doPDFReweight = True

process.pdf3 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf3.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf3.pdfReweightBaseId = 0
process.pdf3.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf3.pdfReweightTargetId = 3
process.pdf3.doPDFReweight = True

process.pdf4 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf4.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf4.pdfReweightBaseId = 0
process.pdf4.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf4.pdfReweightTargetId = 4
process.pdf4.doPDFReweight = True

process.pdf5 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf5.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf5.pdfReweightBaseId = 0
process.pdf5.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf5.pdfReweightTargetId = 5
process.pdf5.doPDFReweight = True

process.pdf6 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf6.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf6.pdfReweightBaseId = 0
process.pdf6.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf6.pdfReweightTargetId = 6
process.pdf6.doPDFReweight = True

process.pdf7 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf7.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf7.pdfReweightBaseId = 0
process.pdf7.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf7.pdfReweightTargetId = 7
process.pdf7.doPDFReweight = True

process.pdf8 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf8.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf8.pdfReweightBaseId = 0
process.pdf8.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf8.pdfReweightTargetId = 8
process.pdf8.doPDFReweight = True

process.pdf9 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf9.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf9.pdfReweightBaseId = 0
process.pdf9.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf9.pdfReweightTargetId = 9
process.pdf9.doPDFReweight = True

process.pdf10 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf10.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf10.pdfReweightBaseId = 0
process.pdf10.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf10.pdfReweightTargetId = 10
process.pdf10.doPDFReweight = True

process.pdf11 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf11.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf11.pdfReweightBaseId = 0
process.pdf11.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf11.pdfReweightTargetId = 11
process.pdf11.doPDFReweight = True

process.pdf12 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf12.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf12.pdfReweightBaseId = 0
process.pdf12.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf12.pdfReweightTargetId = 12
process.pdf12.doPDFReweight = True

process.pdf13 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf13.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf13.pdfReweightBaseId = 0
process.pdf13.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf13.pdfReweightTargetId = 13
process.pdf13.doPDFReweight = True

process.pdf14 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf14.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf14.pdfReweightBaseId = 0
process.pdf14.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf14.pdfReweightTargetId = 14
process.pdf14.doPDFReweight = True

process.pdf15 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf15.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf15.pdfReweightBaseId = 0
process.pdf15.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf15.pdfReweightTargetId = 15
process.pdf15.doPDFReweight = True

process.pdf16 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf16.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf16.pdfReweightBaseId = 0
process.pdf16.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf16.pdfReweightTargetId = 16
process.pdf16.doPDFReweight = True

process.pdf17 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf17.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf17.pdfReweightBaseId = 0
process.pdf17.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf17.pdfReweightTargetId = 17
process.pdf17.doPDFReweight = True

process.pdf18 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf18.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf18.pdfReweightBaseId = 0
process.pdf18.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf18.pdfReweightTargetId = 18
process.pdf18.doPDFReweight = True

process.pdf19 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf19.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf19.pdfReweightBaseId = 0
process.pdf19.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf19.pdfReweightTargetId = 19
process.pdf19.doPDFReweight = True

process.pdf20 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf20.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf20.pdfReweightBaseId = 0
process.pdf20.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf20.pdfReweightTargetId = 20
process.pdf20.doPDFReweight = True

process.pdf21 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf21.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf21.pdfReweightBaseId = 0
process.pdf21.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf21.pdfReweightTargetId = 21
process.pdf21.doPDFReweight = True

process.pdf22 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf22.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf22.pdfReweightBaseId = 0
process.pdf22.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf22.pdfReweightTargetId = 22
process.pdf22.doPDFReweight = True

process.pdf23 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf23.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf23.pdfReweightBaseId = 0
process.pdf23.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf23.pdfReweightTargetId = 23
process.pdf23.doPDFReweight = True

process.pdf24 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf24.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf24.pdfReweightBaseId = 0
process.pdf24.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf24.pdfReweightTargetId = 24
process.pdf24.doPDFReweight = True

process.pdf25 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf25.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf25.pdfReweightBaseId = 0
process.pdf25.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf25.pdfReweightTargetId = 25
process.pdf25.doPDFReweight = True

process.pdf26 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf26.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf26.pdfReweightBaseId = 0
process.pdf26.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf26.pdfReweightTargetId = 26
process.pdf26.doPDFReweight = True

process.pdf27 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf27.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf27.pdfReweightBaseId = 0
process.pdf27.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf27.pdfReweightTargetId = 27
process.pdf27.doPDFReweight = True

process.pdf28 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf28.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf28.pdfReweightBaseId = 0
process.pdf28.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf28.pdfReweightTargetId = 28
process.pdf28.doPDFReweight = True

process.pdf29 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf29.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf29.pdfReweightBaseId = 0
process.pdf29.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf29.pdfReweightTargetId = 29
process.pdf29.doPDFReweight = True

process.pdf30 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf30.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf30.pdfReweightBaseId = 0
process.pdf30.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf30.pdfReweightTargetId = 30
process.pdf30.doPDFReweight = True

process.pdf31 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf31.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf31.pdfReweightBaseId = 0
process.pdf31.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf31.pdfReweightTargetId = 31
process.pdf31.doPDFReweight = True

process.pdf32 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf32.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf32.pdfReweightBaseId = 0
process.pdf32.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf32.pdfReweightTargetId = 32
process.pdf32.doPDFReweight = True

process.pdf33 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf33.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf33.pdfReweightBaseId = 0
process.pdf33.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf33.pdfReweightTargetId = 33
process.pdf33.doPDFReweight = True

process.pdf34 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf34.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf34.pdfReweightBaseId = 0
process.pdf34.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf34.pdfReweightTargetId = 34
process.pdf34.doPDFReweight = True

process.pdf35 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf35.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf35.pdfReweightBaseId = 0
process.pdf35.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf35.pdfReweightTargetId = 35
process.pdf35.doPDFReweight = True

process.pdf36 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf36.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf36.pdfReweightBaseId = 0
process.pdf36.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf36.pdfReweightTargetId = 36
process.pdf36.doPDFReweight = True

process.pdf37 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf37.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf37.pdfReweightBaseId = 0
process.pdf37.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf37.pdfReweightTargetId = 37
process.pdf37.doPDFReweight = True

process.pdf38 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf38.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf38.pdfReweightBaseId = 0
process.pdf38.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf38.pdfReweightTargetId = 38
process.pdf38.doPDFReweight = True

process.pdf39 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf39.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf39.pdfReweightBaseId = 0
process.pdf39.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf39.pdfReweightTargetId = 39
process.pdf39.doPDFReweight = True

process.pdf40 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf40.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf40.pdfReweightBaseId = 0
process.pdf40.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf40.pdfReweightTargetId = 40
process.pdf40.doPDFReweight = True

process.pdf41 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf41.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf41.pdfReweightBaseId = 0
process.pdf41.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf41.pdfReweightTargetId = 41
process.pdf41.doPDFReweight = True

process.pdf42 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf42.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf42.pdfReweightBaseId = 0
process.pdf42.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf42.pdfReweightTargetId = 42
process.pdf42.doPDFReweight = True

process.pdf43 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf43.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf43.pdfReweightBaseId = 0
process.pdf43.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf43.pdfReweightTargetId = 43
process.pdf43.doPDFReweight = True

process.pdf44 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf44.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf44.pdfReweightBaseId = 0
process.pdf44.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf44.pdfReweightTargetId = 44
process.pdf44.doPDFReweight = True

process.pdf45 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf45.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf45.pdfReweightBaseId = 0
process.pdf45.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf45.pdfReweightTargetId = 45
process.pdf45.doPDFReweight = True

process.pdf46 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf46.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf46.pdfReweightBaseId = 0
process.pdf46.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf46.pdfReweightTargetId = 46
process.pdf46.doPDFReweight = True

process.pdf47 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf47.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf47.pdfReweightBaseId = 0
process.pdf47.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf47.pdfReweightTargetId = 47
process.pdf47.doPDFReweight = True

process.pdf48 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf48.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf48.pdfReweightBaseId = 0
process.pdf48.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf48.pdfReweightTargetId = 48
process.pdf48.doPDFReweight = True

process.pdf49 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf49.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf49.pdfReweightBaseId = 0
process.pdf49.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf49.pdfReweightTargetId = 49
process.pdf49.doPDFReweight = True

process.pdf50 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf50.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf50.pdfReweightBaseId = 0
process.pdf50.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf50.pdfReweightTargetId = 50
process.pdf50.doPDFReweight = True

process.pdf51 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf51.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf51.pdfReweightBaseId = 0
process.pdf51.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf51.pdfReweightTargetId = 51
process.pdf51.doPDFReweight = True

process.pdf52 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf52.pdfReweightBaseName = 'CT10.LHgrid'
process.pdf52.pdfReweightBaseId = 0
process.pdf52.pdfReweightTargetName = 'CT10.LHgrid'
process.pdf52.pdfReweightTargetId = 52
process.pdf52.doPDFReweight = True

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
  )

