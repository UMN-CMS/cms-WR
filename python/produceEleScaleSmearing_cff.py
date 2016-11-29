import FWCore.ParameterSet.Config as cms

# make a collection of TuneP muons which pass isHighPt ID
eleScaleSmearing = cms.EDProducer("produceEleScaleSmearing",
                    electrons_src = cms.InputTag("wRminiTreeElectron"),
                    OutputCollectionName1 = cms.string("EleScaleError"),
                    OutputCollectionName2 = cms.string("EleSmearingSigma"),
                    OutputCollectionName3 = cms.string("EleSmearingSigmaPhiUp"),
                    OutputCollectionName4 = cms.string("EleSmearingSigmaPhiDown"),
                    OutputCollectionName5 = cms.string("EleSmearingSigmaRhoUp"),
                    OutputCollectionName6 = cms.string("EleSmearingSigmaRhoDown")
)

