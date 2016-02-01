import FWCore.ParameterSet.Config as cms

MakeTTree_Muons = cms.EDAnalyzer('TTreeMaker',
                            muons_src = cms.InputTag('wRsubleadingMuon'),
                            electrons_src = cms.InputTag('slimmedElectrons'),
                            jets_src = cms.InputTag('wRtightJets'),
                            met_src = cms.InputTag('slimmedMETs'),
                            genparticles_src = cms.InputTag('prunedGenParticles'),
                            genjets_src = cms.InputTag('slimmedGenJets'),
                            primary_vertex_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
                            beamSpot = cms.InputTag('offlineBeamSpot'),
                            isolation_dR = cms.double(0.4),
                            leading_lepton_pt_cut = cms.double(0),
                            lepton_eta_cut = cms.double(5),
                            subleading_lepton_pt_cut = cms.double(0),
                            leading_jet_pt_cut = cms.double(40),
                            jet_eta_cut = cms.double(2.5),
                            subleading_jet_pt_cut = cms.double(40),
                            Mlljj_cut = cms.double(0.0),
                            Mll_cut = cms.double(0.0),
                            muon_mode = cms.bool(True),
                            electron_mode = cms.bool(False),
                            is_mc = cms.bool(True),
                            conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),
                            #
                            # ID decisions (common to all formats)
                            #
                            eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                            eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                            eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                            eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                            eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
                            # An example of configuration for accessing the full cut flow info for
                            # one case is shown below.
                            # The map name for the full info is the same as the map name of the
                            # corresponding simple pass/fail map above, they are distinguished by
                            # the type of the content.
                            eleHEEPIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
                            # This is a fairly verbose mode if switched on, with full cut flow 
                            # diagnostics for each candidate. Use it in a low event count test job.
                            eleIdVerbose = cms.bool(False),
                            #
                            # bTagging
                            #
                            bDiscriminators = cms.vstring(      # list of b-tag discriminators to access
                                'pfTrackCountingHighEffBJetTags',
                                'pfTtrackCountingHighPurBJetTags',
                                'pfJetProbabilityBJetTags',
                                'pfJetBProbabilityBJetTags',
                                'pfSimpleSecondaryVertexHighEffBJetTags',
                                'pfSimpleSecondaryVertexHighPurBJetTags',
                                'pfCombinedSecondaryVertexV2BJetTags',
                                'pfCombinedInclusiveSecondaryVertexV2BJetTags',
                                'pfCombinedMVABJetTags'
                                ),
                            JECUnc_src = cms.InputTag('JECUncertainty'),
                            )
