import FWCore.ParameterSet.Config as cms

MiniTTree = cms.EDAnalyzer('miniTTree',
                            muons_src = cms.InputTag('wRminiTreeMuon'),
                            electrons_src = cms.InputTag('wRminiTreeElectron'),
                            jets_src = cms.InputTag('wRJets'),
                            jec_unc_src = cms.InputTag('wRJECUncert','JECUncertainty'),
                            jetResolution_src = cms.InputTag('jetResolutionSF','JetResolution'),
                            JERsf_src = cms.InputTag('jetResolutionSF','JERsf'),
                            JERsf_up_src = cms.InputTag('jetResolutionSF','JERsfUp'),
                            JERsf_down_src = cms.InputTag('jetResolutionSF','JERsfDown'),
                            muon_IDSF_central_src = cms.InputTag('muonIdIsoSF','MuonSFIdCentral'),                            
                            muon_IsoSF_central_src = cms.InputTag('muonIdIsoSF', 'MuonSFIsoCentral'),
                            muon_IDSF_error_src = cms.InputTag('muonIdIsoSF','MuonSFIdError'),
                            muon_IsoSF_error_src = cms.InputTag('muonIdIsoSF','MuonSFIsoError'),
                            PUWeights_src = cms.InputTag('PUWeights','PileupWeights'),
                           datasetName = cms.InputTag("addStringIdentifier", "datasetIdentifier"),
                            )
