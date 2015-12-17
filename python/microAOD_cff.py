import FWCore.ParameterSet.Config as cms

### \page skim_page Skim 
# (preselection)
# 
# Sequences to produce microAOD are contained in microAOD_cff
#
# it contains the HLT filters and the skim selection sequences
# 
# the output format is defined in microAOD_Output_cff
#
# 
### \par Skim selections
## - leading lepton pt>45 GeV
## - subleading lepton pt>30 GeV
## - jets pt>20 GeV
## - events with 2 leptons and 2 jets (no flavour requirement)
##
## Data and MC are selected using ((OR of all the triggers physics triggers) AND (skim selection)) \b OR (tag&probe)
#
# skims are defined in different cff files:
# - \b skimElectron_cff
# - \b skimMuon_cff
# - \b skimEMu_cff
# - \b skimJets_cff
# and combined in microAOD_cff

""" \addtogroup microAODSeq_Group Sequences to produce microAOD
@{
"""

# set of modules, filters and sequences for skimming
from ExoAnalysis.cmsWR.skim_cff import *

# hlt trigger paths used in the analysis
from ExoAnalysis.cmsWR.hltFilters_cff import *


############################################################

# # Pruning the collections with only preselected electrons
# from PhysicsTools.PatAlgos.slimming.slimmedElectrons_cfi import *
# slimmedWrSelectedElectrons = slimmedElectrons.clone()
# #slimmedWrSelectedElectrons.src=cms.InputTag('wRpreSelectedElectrons')
# slimmedWrSelectedElectrons.src=cms.InputTag('wRleadingElectron')
# slimmedWrSelectedElectrons.linkToPackedPFCandidates = cms.bool(False)
# slimmedWrSelectedElectrons.reducedBarrelRecHitCollection = cms.InputTag("reducedEgamma","reducedEBRecHits")
# slimmedWrSelectedElectrons.reducedEndcapRecHitCollection = cms.InputTag("reducedEgamma","reducedEERecHits")

#wRslimmingElectronSeq = cms.Sequence(slimmedWrSelectedElectrons)

#microAODslimmingSeq = cms.Sequence(wRslimmingElectronSeq)
"""
@}
"""
