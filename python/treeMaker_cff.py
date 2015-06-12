import FWCore.ParameterSet.Config as cms


# Sequences to produce microAOD
""" \addtogroup microAODSeq_Group Sequences to produce microAOD
@{
"""

from ExoAnalysis.cmsWR.skimMuon_cff import *
from ExoAnalysis.cmsWR.skimElectron_cff import *

TFileService = cms.Service('TFileService', fileName = cms.string('made_ttree.root'))

MakeTTree = cms.EDAnalyzer('TTreeMaker',
                            muons_src = cms.InputTag('wRsubleadingMuons'),
                            electrons_src = cms.InputTag('wRsubleadingElectron'),
                            jets_src = cms.InputTag('slimmedJets'),
                            genparticles_src = cms.InputTag('prunedGenParticles'),
                            genjets_src = cms.InputTag('slimmedGenJets'),
                            primary_vertex_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
                            isolation_dR = cms.double(0.0),
                            leading_lepton_pt_cut = cms.double(0),
                            lepton_eta_cut = cms.double(5),
                            subleading_lepton_pt_cut = cms.double(0),
                            leading_jet_pt_cut = cms.double(0),
                            jet_eta_cut = cms.double(5),
                            subleading_jet_pt_cut = cms.double(0),
                            Mlljj_cut = cms.double(0.0),
                            Mll_cut = cms.double(0.0),
                            muon_mode = cms.bool(True)
                            )

"""
@}
"""
