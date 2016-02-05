import FWCore.ParameterSet.Config as cms

MiniTTree = cms.EDAnalyzer('miniTTree',
                            muons_src = cms.InputTag('wRminiTreeMuon'),
                            electrons_src = cms.InputTag('wRminiTreeElectron'),
                            jets_src = cms.InputTag('wRJets'),
                            jec_unc_src = cms.InputTag('wRJECUncert','JECUncertainty'),
                            muon_IDSF_central_src = cms.InputTag('muonIdIsoSF','MuonSFIdCentral'),
                            muon_IsoSF_central_src = cms.InputTag('muonIdIsoSF', 'MuonSFIsoCentral'),
                            muon_IDSF_error_src = cms.InputTag('muonIdIsoSF','MuonSFIdError'),
                            muon_IsoSF_error_src = cms.InputTag('muonIdIsoSF','MuonSFIsoError'),
                            PUWeights_src = cms.InputTag('PUWeights','PileupWeights'),
                            )
