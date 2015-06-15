#use this file to define the PoolSource input and other modules which will be used
#repeatedly in skim_cfg.py
import FWCore.ParameterSet.Config as cms

source = cms.Source("PoolSource",
		fileNames = cms.untracked.vstring('file:noFile.root'),
		secondaryFileNames = cms.untracked.vstring()
		)

