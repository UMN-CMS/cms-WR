import FWCore.ParameterSet.Config as cms

# make a collection of TuneP muons which pass isHighPt ID
ScaleCorrectedMuonsProd = cms.EDProducer("produceScaleCorrectedMuons",
		src = cms.InputTag("slimmedMuons"),
                OutputCollectionName = cms.string("ScaleCorrectedMuons")
		)

ScaleCorrectedMuonsSequence = cms.Sequence(ScaleCorrectedMuonsProd)
