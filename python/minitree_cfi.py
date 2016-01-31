import FWCore.ParameterSet.Config as cms

MiniTTree = cms.EDAnalyzer('miniTTree',
                            muons_src = cms.InputTag('wRsubleadingMuon'),
                            electrons_src = cms.InputTag('wRminiTreeElectron'),
                            jets_src = cms.InputTag('wRJets'),
                            jec_unc_src = cms.InputTag('wRJECUncert','JECUncertainty'),
                            muon_IDSF_central_src = cms.InputTag('MuonSFIdCentral','muonIdIsoSF'),
                            muon_IsoSF_central_src = cms.InputTag('MuonSFIsoCentral','muonIdIsoSF'),
                            muon_IDSF_error_src = cms.InputTag('MuonSFIdError','muonIdIsoSF'),
                            muon_IsoSF_error_src = cms.InputTag('MuonSFIsoError','muonIdIsoSF'),
                            )
