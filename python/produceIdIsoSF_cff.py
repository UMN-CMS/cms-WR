import FWCore.ParameterSet.Config as cms

import pickle
f = open('/afs/cern.ch/work/r/rchatter/CMSSW_7_4_15_patch1/src/ExoAnalysis/cmsWR/python/MuonHighPt_Z_RunCD_Reco74X_Dec17.pkl', 'r')
results = pickle.load(f)
results.keys()
results["HighPtID_EtaBins_Pt53"].keys()
ii = 0
SF_ID_C = []
SF_ID_E = []
SF_ISO_C = []
SF_ISO_E = []

for i in xrange(0,13):
    SF_ID_C.append(0)
    SF_ID_E.append(0)
    SF_ISO_C.append(0)
    SF_ISO_E.append(0)

for key, result in sorted(results["HighPtID_EtaBins_Pt53"]["pTtuneP_ratio"].iteritems()) : 
    SF_ID_C[ii] = result["value"]
    SF_ID_E[ii] = result["error"]
    ii += 1

ii=0
for key, result in sorted(results["tkRelIsoID_EtaBins_Pt53"]["eta_ratio"].iteritems()) :
    SF_ISO_C[ii] = result["value"]
    SF_ISO_E[ii] = result["error"]
    ii += 1

    

# make a collection of TuneP muons which pass isHighPt ID
MuonIdIsoSFProd = cms.EDProducer("produceLepIdIsoScaleFactors",
		src = cms.InputTag("ScaleCorrectedMuonsProd","ScaleCorrectedMuons"),
                OutputCollectionName1 = cms.string("MuonSFIdCentral"),
                OutputCollectionName2 = cms.string("MuonSFIdError"),
                OutputCollectionName3 = cms.string("MuonSFIsoCentral"),
                OutputCollectionName4 = cms.string("MuonSFIsoError"),
                Scale_Factor_ID_Central = cms.vdouble(SF_ID_C),
                Scale_Factor_ID_Error = cms.vdouble(SF_ID_E),
                Scale_Factor_ISO_Central = cms.vdouble(SF_ISO_C),
                Scale_Factor_ISO_Error = cms.vdouble(SF_ISO_E)
)

MuonIdIsoSFProdSequence = cms.Sequence(MuonIdIsoSFProd)
