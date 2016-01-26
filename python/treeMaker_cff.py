import FWCore.ParameterSet.Config as cms


# TTree producer
""" \addtogroup TTree producer
@{
"""

from ExoAnalysis.cmsWR.treeMaker_cfi import *

MakeTTree_Electrons = MakeTTree_Muons.clone(muon_mode = cms.bool(False),
                                            electron_mode = cms.bool(True))
MakeTTree_EleMuon = MakeTTree_Muons.clone(muon_mode = cms.bool(True),
                                            electron_mode = cms.bool(True))
"""
@}
"""
