import FWCore.ParameterSet.Config as cms

# Pileup root files should contain a TH1F in its rootdir named "pileup"
PUWeights = cms.EDProducer("PileupWeight",
		outputCollectionName = cms.string("PileupWeights"),
		PileupMCFilename = cms.string("data/MCPileup.root"),
		PileupDataFilename = cms.string("data/DataPileup.root"),
)

