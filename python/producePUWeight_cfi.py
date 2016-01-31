import FWCore.ParameterSet.Config as cms

# Pileup root files should contain a TH1F in its rootdir named "pileup"
PUWeights = cms.EDProducer("producePileupWeight",
		outputCollectionName = cms.string("PileupWeights"),
		PileupMCFilename = cms.string("MCPileup.root"),
		PileupDataFilename = cms.string("DataPileup.root"),
)

PUWeightsSequence = cms.Sequence(PUWeights)
