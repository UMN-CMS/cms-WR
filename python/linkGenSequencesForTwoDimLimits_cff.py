import FWCore.ParameterSet.Config as cms

from ExoAnalysis.cmsWR.genEleChnlCutSeqForTwoDimLimits_cff import *
from ExoAnalysis.cmsWR.genMuMuChnlCutSeqForTwoDimLimits_cff import *

#setup the evt count producers which will be run before any cuts and after all cuts
countEvtsBeforeCuts = cms.EDProducer("produceEventCount",
		outputCollectionName = cms.string("evtCountBeforeAllCuts")
		)

countEvtsAfterCuts = cms.EDProducer("produceEventCount",
		outputCollectionName = cms.string("evtCountAfterAllCuts")
		)

#setup the analyzer which stores the number of evts before and after cuts, the cut efficiency,
#and wr and nu mass in a TTree
cutEfficiencyAndWrAndNuMasses = cms.EDAnalyzer("determineCutEfficiency",
		massOfWR = cms.double(2000),
		massOfNu = cms.double(1000),
		treeName = cms.string("cutEffTree"),
		evtCountBeforeCutsCollection = cms.InputTag("countEvtsBeforeCuts","evtCountBeforeAllCuts"),
		evtCountAfterCutsCollection = cms.InputTag("countEvtsAfterCuts","evtCountAfterAllCuts")
		)

#eeChnlCutSeq applies all of the offline cuts using gen electrons and quarks from the genParticles collection
#the specific cut values can be changed within genEleChnlCutSeqForTwoDimLimits_cff.py
eeChnlCutSeq = cms.Sequence(
		countEvtsBeforeCuts
		*eeChnlPtEtaFilteredMatchedGenParticleSeq
		*eeChnlMatchedDiLeptonCandidateSeq
		*eeChnlMatchedQuarkLeptonDrSepSeq
		*eeChnlMatchedFourObjectCandidateSeq
		*countEvtsAfterCuts
		)

#mumuChnlCutSeq applies all of the offline cuts using gen muons and quarks from the genParticles collection
#the specific cut values can be changed within genMuMuChnlCutSeqForTwoDimLimits_cff.py
mumuChnlCutSeq = cms.Sequence(
		countEvtsBeforeCuts
		*mumuChnlPtEtaFilteredMatchedGenParticleSeq
		*mumuChnlMatchedDiLeptonCandidateSeq
		*mumuChnlMatchedQuarkLeptonDrSepSeq
		*mumuChnlMatchedFourObjectCandidateSeq
		*countEvtsAfterCuts
		)

