import FWCore.ParameterSet.Config as cms

MiniTTree = cms.EDAnalyzer('miniTTree',
                            muons_src = cms.InputTag('wRsubleadingMuon'),
                            electrons_src = cms.InputTag('wRminiTreeElectron'),
                            jets_src = cms.InputTag('wRJets'),
                            jec_unc_src = cms.InputTag('wRJECUncert'),
                            muon_IDSF_central_src = cms.InputTag('muonIdIsoSF','MuonSFIdCentral'),
                           muon_IsoSF_central_src = cms.InputTag('muonIdIsoSF', 'MuonSFIsoCentral'),
#                OutputCollectionName2 = cms.string("MuonSFIdError"),
#                OutputCollectionName4 = cms.string("MuonSFIsoError"),

                            )
