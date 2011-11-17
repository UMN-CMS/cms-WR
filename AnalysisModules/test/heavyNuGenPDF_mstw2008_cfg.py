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
process.pdf0.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf0.pdfReweightBaseId = 0
process.pdf0.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf0.pdfReweightTargetId = 0
process.pdf0.doPDFReweight = True

process.pdf1 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf1.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf1.pdfReweightBaseId = 0
process.pdf1.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf1.pdfReweightTargetId = 1
process.pdf1.doPDFReweight = True

process.pdf2 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf2.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf2.pdfReweightBaseId = 0
process.pdf2.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf2.pdfReweightTargetId = 2
process.pdf2.doPDFReweight = True

process.pdf3 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf3.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf3.pdfReweightBaseId = 0
process.pdf3.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf3.pdfReweightTargetId = 3
process.pdf3.doPDFReweight = True

process.pdf4 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf4.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf4.pdfReweightBaseId = 0
process.pdf4.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf4.pdfReweightTargetId = 4
process.pdf4.doPDFReweight = True

process.pdf5 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf5.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf5.pdfReweightBaseId = 0
process.pdf5.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf5.pdfReweightTargetId = 5
process.pdf5.doPDFReweight = True

process.pdf6 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf6.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf6.pdfReweightBaseId = 0
process.pdf6.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf6.pdfReweightTargetId = 6
process.pdf6.doPDFReweight = True

process.pdf7 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf7.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf7.pdfReweightBaseId = 0
process.pdf7.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf7.pdfReweightTargetId = 7
process.pdf7.doPDFReweight = True

process.pdf8 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf8.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf8.pdfReweightBaseId = 0
process.pdf8.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf8.pdfReweightTargetId = 8
process.pdf8.doPDFReweight = True

process.pdf9 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf9.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf9.pdfReweightBaseId = 0
process.pdf9.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf9.pdfReweightTargetId = 9
process.pdf9.doPDFReweight = True

process.pdf10 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf10.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf10.pdfReweightBaseId = 0
process.pdf10.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf10.pdfReweightTargetId = 10
process.pdf10.doPDFReweight = True

process.pdf11 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf11.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf11.pdfReweightBaseId = 0
process.pdf11.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf11.pdfReweightTargetId = 11
process.pdf11.doPDFReweight = True

process.pdf12 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf12.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf12.pdfReweightBaseId = 0
process.pdf12.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf12.pdfReweightTargetId = 12
process.pdf12.doPDFReweight = True

process.pdf13 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf13.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf13.pdfReweightBaseId = 0
process.pdf13.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf13.pdfReweightTargetId = 13
process.pdf13.doPDFReweight = True

process.pdf14 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf14.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf14.pdfReweightBaseId = 0
process.pdf14.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf14.pdfReweightTargetId = 14
process.pdf14.doPDFReweight = True

process.pdf15 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf15.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf15.pdfReweightBaseId = 0
process.pdf15.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf15.pdfReweightTargetId = 15
process.pdf15.doPDFReweight = True

process.pdf16 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf16.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf16.pdfReweightBaseId = 0
process.pdf16.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf16.pdfReweightTargetId = 16
process.pdf16.doPDFReweight = True

process.pdf17 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf17.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf17.pdfReweightBaseId = 0
process.pdf17.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf17.pdfReweightTargetId = 17
process.pdf17.doPDFReweight = True

process.pdf18 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf18.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf18.pdfReweightBaseId = 0
process.pdf18.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf18.pdfReweightTargetId = 18
process.pdf18.doPDFReweight = True

process.pdf19 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf19.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf19.pdfReweightBaseId = 0
process.pdf19.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf19.pdfReweightTargetId = 19
process.pdf19.doPDFReweight = True

process.pdf20 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf20.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf20.pdfReweightBaseId = 0
process.pdf20.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf20.pdfReweightTargetId = 20
process.pdf20.doPDFReweight = True

process.pdf21 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf21.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf21.pdfReweightBaseId = 0
process.pdf21.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf21.pdfReweightTargetId = 21
process.pdf21.doPDFReweight = True

process.pdf22 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf22.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf22.pdfReweightBaseId = 0
process.pdf22.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf22.pdfReweightTargetId = 22
process.pdf22.doPDFReweight = True

process.pdf23 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf23.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf23.pdfReweightBaseId = 0
process.pdf23.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf23.pdfReweightTargetId = 23
process.pdf23.doPDFReweight = True

process.pdf24 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf24.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf24.pdfReweightBaseId = 0
process.pdf24.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf24.pdfReweightTargetId = 24
process.pdf24.doPDFReweight = True

process.pdf25 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf25.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf25.pdfReweightBaseId = 0
process.pdf25.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf25.pdfReweightTargetId = 25
process.pdf25.doPDFReweight = True

process.pdf26 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf26.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf26.pdfReweightBaseId = 0
process.pdf26.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf26.pdfReweightTargetId = 26
process.pdf26.doPDFReweight = True

process.pdf27 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf27.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf27.pdfReweightBaseId = 0
process.pdf27.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf27.pdfReweightTargetId = 27
process.pdf27.doPDFReweight = True

process.pdf28 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf28.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf28.pdfReweightBaseId = 0
process.pdf28.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf28.pdfReweightTargetId = 28
process.pdf28.doPDFReweight = True

process.pdf29 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf29.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf29.pdfReweightBaseId = 0
process.pdf29.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf29.pdfReweightTargetId = 29
process.pdf29.doPDFReweight = True

process.pdf30 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf30.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf30.pdfReweightBaseId = 0
process.pdf30.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf30.pdfReweightTargetId = 30
process.pdf30.doPDFReweight = True

process.pdf31 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf31.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf31.pdfReweightBaseId = 0
process.pdf31.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf31.pdfReweightTargetId = 31
process.pdf31.doPDFReweight = True

process.pdf32 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf32.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf32.pdfReweightBaseId = 0
process.pdf32.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf32.pdfReweightTargetId = 32
process.pdf32.doPDFReweight = True

process.pdf33 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf33.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf33.pdfReweightBaseId = 0
process.pdf33.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf33.pdfReweightTargetId = 33
process.pdf33.doPDFReweight = True

process.pdf34 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf34.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf34.pdfReweightBaseId = 0
process.pdf34.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf34.pdfReweightTargetId = 34
process.pdf34.doPDFReweight = True

process.pdf35 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf35.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf35.pdfReweightBaseId = 0
process.pdf35.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf35.pdfReweightTargetId = 35
process.pdf35.doPDFReweight = True

process.pdf36 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf36.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf36.pdfReweightBaseId = 0
process.pdf36.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf36.pdfReweightTargetId = 36
process.pdf36.doPDFReweight = True

process.pdf37 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf37.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf37.pdfReweightBaseId = 0
process.pdf37.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf37.pdfReweightTargetId = 37
process.pdf37.doPDFReweight = True

process.pdf38 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf38.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf38.pdfReweightBaseId = 0
process.pdf38.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf38.pdfReweightTargetId = 38
process.pdf38.doPDFReweight = True

process.pdf39 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf39.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf39.pdfReweightBaseId = 0
process.pdf39.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf39.pdfReweightTargetId = 39
process.pdf39.doPDFReweight = True

process.pdf40 =HeavyNu.AnalysisModules.heavynugenlevel_cfi.hNuGen.clone()
process.pdf40.pdfReweightBaseName = 'MSTW2008lo68cl.LHgrid'
process.pdf40.pdfReweightBaseId = 0
process.pdf40.pdfReweightTargetName = 'MSTW2008lo68cl.LHgrid'
process.pdf40.pdfReweightTargetId = 40
process.pdf40.doPDFReweight = True

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
  )

