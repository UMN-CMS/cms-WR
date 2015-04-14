import FWCore.ParameterSet.Config as cms


# Sequences to produce microAOD
""" \addtogroup microAODSeq_Group Sequences to produce microAOD
@{
"""

from ExoAnalysis.cmsWR.skimMuon_cff import *
from ExoAnalysis.cmsWR.skimElectron_cff import *


# Pruning the collections with only preselected electrons
from PhysicsTools.PatAlgos.slimming.slimmedElectrons_cfi import *
slimmedWrSelectedElectrons = slimmedElectrons.clone()
#slimmedWrSelectedElectrons.src=cms.InputTag('wRpreSelectedElectrons')
slimmedWrSelectedElectrons.src=cms.InputTag('wRleadingElectron')
slimmedWrSelectedElectrons.linkToPackedPFCandidates = cms.bool(False)
slimmedWrSelectedElectrons.reducedBarrelRecHitCollection = cms.InputTag("reducedEgamma","reducedEBRecHits")
slimmedWrSelectedElectrons.reducedEndcapRecHitCollection = cms.InputTag("reducedEgamma","reducedEERecHits")

wRslimmingElectronSeq = cms.Sequence(slimmedWrSelectedElectrons)

microAODslimmingSeq = cms.Sequence(wRslimmingElectronSeq)
"""
@}
"""
