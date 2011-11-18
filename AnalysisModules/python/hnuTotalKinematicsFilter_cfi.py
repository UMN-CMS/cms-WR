import FWCore.ParameterSet.Config as cms

hnuTotalKinematicsFilter = cms.EDFilter('HnuTotalKinematicsFilter',
  src             = cms.InputTag("genParticles"),
  tolerance       = cms.double(0.5),
  verbose         = cms.untracked.bool(False)                                   
)
