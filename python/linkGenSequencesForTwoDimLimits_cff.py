import FWCore.ParameterSet.Config as cms

from ExoAnalysis.cmsWR.genEleChnlCutSeqForTwoDimLimits_cff import *
from ExoAnalysis.cmsWR.genMuMuChnlCutSeqForTwoDimLimits_cff import *


#eeChnlCutSeq applies all of the offline cuts using gen electrons and quarks from the genParticles collection
#the specific cut values can be changed within genEleChnlCutSeqForTwoDimLimits_cff.py
eeChnlCutSeq = cms.Sequence(
		eeChnlPtEtaFilteredMatchedGenParticleSeq
		*eeChnlMatchedDiLeptonCandidateSeq
		*eeChnlMatchedQuarkLeptonDrSepSeq
		*eeChnlMatchedFourObjectCandidateSeq
		)

#mumuChnlCutSeq applies all of the offline cuts using gen muons and quarks from the genParticles collection
#the specific cut values can be changed within genMuMuChnlCutSeqForTwoDimLimits_cff.py
mumuChnlCutSeq = cms.Sequence(
		mumuChnlPtEtaFilteredMatchedGenParticleSeq
		*mumuChnlMatchedDiLeptonCandidateSeq
		*mumuChnlMatchedQuarkLeptonDrSepSeq
		*mumuChnlMatchedFourObjectCandidateSeq
		)

