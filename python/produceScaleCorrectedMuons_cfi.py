import FWCore.ParameterSet.Config as cms

# make a collection of TuneP muons which pass isHighPt ID
ScaleCorrectedMuonsProd = cms.EDProducer("produceScaleCorrectedMuons",
                                         src = cms.InputTag("tunePIDIsoMuons"),
                                         OutputCollectionName = cms.string("ScaleCorrectedtunePIDIsoMuons"),
                                         initFile = cms.string("ExoAnalysis/cmsWR/data/RoccoR_13tev.txt"),
)
