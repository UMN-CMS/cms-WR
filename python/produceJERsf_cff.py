import FWCore.ParameterSet.Config as cms

# make a collection of TuneP muons which pass isHighPt ID
jetResolutionSF = cms.EDProducer("produceJetResolution",
                    jets_src = cms.InputTag("wRJets"),
                    rho_src = cms.InputTag("fixedGridRhoAll"),
                    OutputCollectionName1 = cms.string("JetResolution"),
                    OutputCollectionName2 = cms.string("JERsf"),
                    OutputCollectionName3 = cms.string("JERsfUp"),
                    OutputCollectionName4 = cms.string("JERsfDown")
)

