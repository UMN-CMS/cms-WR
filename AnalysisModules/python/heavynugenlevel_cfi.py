import FWCore.ParameterSet.Config as cms

hNuGen = cms.EDAnalyzer (
    "HeavyNuGenLevel",
    doPDFReweight = cms.untracked.bool(False),
#    pdfReweightBaseName = cms.untracked.string('cteq6ll.LHpdf'),
    pdfReweightBaseName = cms.untracked.string('CT10.LHgrid'),
    pdfReweightBaseId = cms.untracked.int32(0),
    pdfReweightTargetName = cms.untracked.string('CT10'),
    pdfReweightTargetId = cms.untracked.int32(0),
    pdfReweightAddZMass = cms.untracked.bool(False),
    minMu1pt     = cms.double(60.),
    minMu2pt     = cms.double(40.),
    minJetPt     = cms.double(40),
    maxMuAbsEta  = cms.double(2.4),
    maxJetAbsEta = cms.double(2.5),
    minMuMuMass  = cms.double(200),
    min4objMass  = cms.double(600),
    minMuonJetdR = cms.double(0.5),
    )
