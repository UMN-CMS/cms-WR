import FWCore.ParameterSet.Config as cms

MiniTTree = cms.EDAnalyzer('miniTTree',
                            muons_src = cms.InputTag('wRsubleadingMuon'),
                            electrons_src = cms.InputTag('wRsubleadingElectron'),
                            jets_src = cms.InputTag('wRtightJets'),
                            is_mc = cms.bool(False),
                            )
