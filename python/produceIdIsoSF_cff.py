import FWCore.ParameterSet.Config as cms

import pickle
f = open('/afs/cern.ch/work/r/rchatter/CMSSW_7_4_15_patch1/src/ExoAnalysis/cmsWR/python/MuonHighPt_Z_RunCD_Reco74X_Dec17.pkl', 'r')
results = pickle.load(f)
results.keys()
results["HighPtID_EtaBins_Pt53"].keys()
ii = 0
SF_C = []
SF_E = []

for i in xrange(0,13):
    SF_C.append(0)
    SF_E.append(0)

for key, result in sorted(results["HighPtID_EtaBins_Pt53"]["pTtuneP_ratio"].iteritems()) : 
    SF_C[ii] = result["value"]
    SF_E[ii] = result["error"]
    ii += 1

    

# make a collection of TuneP muons which pass isHighPt ID
MuonIdIsoSFProd = cms.EDProducer("produceLepIdIsoScaleFactors",
		src = cms.InputTag("slimmedMuons"),
                OutputCollectionName1 = cms.string("SFIdIsoCentral"),
                OutputCollectionName2 = cms.string("SFIdIsoError"),
                Scale_Factor_Central = cms.vdouble(SF_C),
                Scale_Factor_Error = cms.vdouble(SF_E)
)

MuonIdIsoSFProdSequence = cms.Sequence(MuonIdIsoSFProd)
