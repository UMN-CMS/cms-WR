import FWCore.ParameterSet.Config as cms
from ExoAnalysis.cmsWR.pileupWeight_cfi import *

PUWeightsSequence = cms.Sequence(PUWeights)
