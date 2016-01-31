import FWCore.ParameterSet.Config as cms

# make a collection of TuneP muons which pass isHighPt ID
ScaleCorrectedMuonsProd = cms.EDProducer("produceScaleCorrectedMuons",
		src = cms.InputTag("tunePIDIsoMuons"),
                OutputCollectionName = cms.string("ScaleCorrectedtunePIDIsoMuons")
		)

ScaleCorrectedMuonsSequence = cms.Sequence(ScaleCorrectedMuonsProd)
