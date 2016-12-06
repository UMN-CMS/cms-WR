import FWCore.ParameterSet.Config as cms

checkWeights = cms.EDAnalyzer('factAndRenormWeightAnalyzer',
		treeName = cms.string("weightsTree")
		)
