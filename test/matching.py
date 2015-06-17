import os,glob,sys
import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('root://xrootd.unl.edu//store/mc/RunIISpring15DR74/WRToNuMuToMuMuJJ_MW-2600_MNu-1300_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/102CF7B3-1F08-E511-93BB-00074305CF52.root'),
                            )
        
outfile = 'matching.root'

process.TFileService = cms.Service('TFileService', fileName = cms.string(outfile))

#from ExoAnalysis.cmsWR.skimMuon_cff import *
#from ExoAnalysis.cmsWR.skimElectron_cff import *
process.load('ExoAnalysis.cmsWR.skimMuon_cff')


process.matching = cms.EDAnalyzer('Matching',
                            muon_src = cms.InputTag('wRsubleadingMuon'),
                            vertex_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
                            electron_src = cms.InputTag('wRsubleadingElectron'),
                            jet_src = cms.InputTag('slimmedJets'),
                            gen_src = cms.InputTag('prunedGenParticles'),
                            gen_jet_src = cms.InputTag('slimmedGenJets')
                            )


process.p = cms.Path(process.wRtunePMuons + process.wRsubleadingMuon + process.matching)
