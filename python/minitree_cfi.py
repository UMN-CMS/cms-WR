import FWCore.ParameterSet.Config as cms

MiniTTree = cms.EDAnalyzer('miniTTree',
                            muons_src = cms.InputTag('wRsubleadingMuon'),
                            electrons_src = cms.InputTag('wRminiTreeElectron'),
                            jets_src = cms.InputTag('wRJets'),
                            jec_unc_src = cms.string('wRJECUncert'),
                            muon_IDSF_central_src = cms.string(''),
                            muon_IsoSF_central_src = cms.string(''),
                            )
