import FWCore.ParameterSet.Config as cms

MiniTTree = cms.EDAnalyzer('miniTTree',
                            muons_src = cms.InputTag('wRsubleadingMuon'),
                            electrons_src = cms.InputTag('wRminiTreeElectron'),
                            jets_src = cms.InputTag('wRJets'),
                            jec_unc_src = cms.string('wRECUncert'),                            
                            )
