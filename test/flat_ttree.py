import sys
import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('file:EXO-Phys14DR-00009.root'),
                            )

process.TFileService = cms.Service('TFileService', fileName = cms.string('flat_ttree.root'))

process.wRtunePMuons = cms.EDProducer("TunePMuonProducer",
                             src = cms.InputTag("slimmedMuons")
                             )

process.ft_nocuts = cms.EDAnalyzer('FlatTTree',
                            muons_src = cms.InputTag('wRtunePMuons'),
                            jets_src = cms.InputTag('slimmedJets'),
                            genparticles_src = cms.InputTag('prunedGenParticles'),
                            genjets_src = cms.InputTag('slimmedGenJets'),
                            primary_vertex_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
                            Mlljj_cut = cms.double(0.0),
                            Mll_cut = cms.double(0.0),
                            matching_dR = cms.double(0.4),
                            matching = cms.bool(True),
                            )

process.p = cms.Path(process.wRtunePMuons + process.ft_nocuts)
